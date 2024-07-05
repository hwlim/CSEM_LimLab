#!/usr/bin/env python3

'''
    Written by:
        Christopher Ahn

    Script used to pick one read pair (among multimappers) based on the posterior probability noted as "ZW:f:n" in the tags, where n is the probability
    Input: bam file
        - needs "ZW:f:n" information in the tags section where n is the probability
        - bam file needs to be sorted by read name
        - no bam index required
    Output: bam file containing only one read pair each
'''

#import required modules and check if user has pysam installed
try:
    import pysam
except ImportError:
    sys.exit("Error: This script requires the module \"pysam\" to run.")
import argparse
import sys
import random

#parse command line arguments
options = argparse.ArgumentParser(description="Unify BAM file from CSEM by random selection following CSEM weights. Picks a single pair of reads among multimappers from a read name-sorted bam file based on a posterior probability. Posterior probability should be denoted in the metadata as \"ZW:f:n\", where n is the probability. Unique reads will have a posterior probabilty of 1, thus will be output as is.", usage="ngs.pickFromMultimappers.py (options) [bam]")
options.add_argument('-o', '--outFile', default='output.bam',
                        help='Output file name; default = output.bam')
options.add_argument('bam',
                        help='Required; Read sorted bam file containing \"ZW:f:n\" in the metadata; n is the posterior probability.')
args = options.parse_args()

## takes list of all reads of the same read name and a list of probabilities for each pair
## return a pair of reads in the form of a list (always contains 2 items)
def pick_pair(reads, probs):

    ## convert zero probability to 1e-99
    probs = [1e-99 if prob == 0 else prob for prob in probs]

    finalPair = random.choices(reads, probs)[0]
    
    return finalPair

## get posterior probability; exit if ZW tag doesn't exist
def get_prob(read):
    
    ## Iterate thru tags backwards, since ZW will probably be towards the end
    for tag, value in reversed(read.tags):
        if tag.startswith("ZW"):
            return float(value)
        
    ## If ZW not found, exit
    print('"ZW" tag is not found in the metadata section of the bam file. Exiting...')
    exit()

## main function
def main():

    ## initialize variables
    savedReadName = ""

    ## list containing a pair; always length of 2
    readPairList = []

    ## list containing readPairLists; length n, where n is from NH:i:n
    allReadsList = []

    ## list containing posterior probabilities; length n, where n is from NH:i:n
    readProbList = []

    ## used to distinguish first and second pair
    i = 0

    ## open input file
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        
        ## open output file
        with pysam.AlignmentFile(args.outFile, "wb", header=bam.header) as bamOutputFile:
            
            ## go through each line of input file
            for read in bam:
                
                ## get readName
                readName = read.query_name

                ## get posterior probability
                posterior = get_prob(read)
                
                ## handle first line of input bam
                if i == 0:
                    savedReadName = readName
                    readPairList.append(read)
                    readProbList.append(posterior)
                    i += 1
                    continue

                ## if this read name is not seen for the first time
                if readName == savedReadName:

                    ## if this is the first of the pair
                    if i % 2 == 0:

                        allReadsList.append(readPairList)
                        readProbList.append(posterior)
                        readPairList = [read]
                    
                    else:
                        readPairList.append(read)

                ## if this read is seen for the first time
                else:

                    ## finish constructing list of reads
                    allReadsList.append(readPairList)

                    ## if unique alignment, don't do random.choices()
                    if len(readProbList) == 1:
                        finalPair = allReadsList[0]
                    
                    ## else do random.choices()
                    else:
                        finalPair = pick_pair(allReadsList, readProbList)
                    
                    ## write output to bam file
                    bamOutputFile.write(finalPair[0])
                    bamOutputFile.write(finalPair[1])
                    
                    ## reset variables
                    readPairList = [read]
                    allReadsList = []
                    readProbList = [posterior]
                    savedReadName = readName
                
                i += 1

            ## take care of last pair of input bam file
            allReadsList.append(readPairList)

            finalPair = pick_pair(allReadsList, readProbList)
            
            ## write output to bam file
            bamOutputFile.write(finalPair[0])
            bamOutputFile.write(finalPair[1])

if __name__ == "__main__":
    main()