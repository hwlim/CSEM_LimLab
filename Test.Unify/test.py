#!/usr/bin/env python3

'''
    Written by:
        Christopher Ahn

    Script used to test and check frequencies of unifyCSEM and output a scatterplot.
    Input: test bam file
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
import os
import re
import sys
import random
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt


#parse command line arguments
options = argparse.ArgumentParser(description="Test ngs.unifyCSEM.py by sampling reads in a test bam file n times and showing the results as a scatterplot.", usage="test.py (options) [bam]")
options.add_argument('-o', '--outPrefix', default='~/scatterplot',
                        help='Output file name for scatterplot; default = output.[png,pdf]')
options.add_argument('-n', '--repeats', default=10000,
                        help='Number of times to repeat the simulation; default = 10000')
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

## draw scatterplot after sampling
def main():

    n = int(args.repeats)
    bamIn = args.bam
    outPrefix = args.outPrefix

    ## sample test
    samplingList = []

    ## initialize variables
    savedReadName = ""

    ## list containing a pair; always length of 2
    readPairList = []

    ## list containing readPairLists; length n, where n is from NH:i:n
    allReadsList = []

    ## list containing posterior probabilities; length n, where n is from NH:i:n
    readProbList = []

    for j in range(n):
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
        with pysam.AlignmentFile(bamIn, "rb") as bam:
        
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

                    ## print
                    item = f'{finalPair[0].query_name}-{bam.get_reference_name(finalPair[0].reference_id)}-{finalPair[0].reference_start}-{finalPair[0].tags[-1][-1]}'
                    samplingList.append(item)
                    
                    ## reset variables
                    readPairList = [read]
                    allReadsList = []
                    readProbList = [posterior]
                    savedReadName = readName
                
                i += 1
        
            ## take care of last pair of input bam file
            allReadsList.append(readPairList)
        
            finalPair = pick_pair(allReadsList, readProbList)
            
            ## print
            item = f'{finalPair[0].query_name}-{bam.get_reference_name(finalPair[0].reference_id)}-{finalPair[0].reference_start}-{finalPair[0].tags[-1][-1]}'
            samplingList.append(item)


    # Create a pandas Series
    series_data = pd.Series(samplingList)

    # Generate the frequency table
    frequency_table = series_data.value_counts()

    df = frequency_table.reset_index()
    split_columns = df['index'].str.split('-', expand=True)

    # Optionally, rename the new columns if needed
    split_columns.columns = ['Col1', 'Col2', 'Col3', 'Col4']
    df = pd.concat([split_columns, df.drop(columns=['index'])], axis=1)

    # Extract the relevant columns for the scatter plot
    x_values = df['Col4'].astype(float)  # Second last column
    y_values = df['count']  # Last column

    # Plotting the scatterplot
    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, y_values, c='blue', alpha=0.5)

    # Setting the x-axis range and ticks
    plt.xlim(0, 1)
    plt.xticks(rotation=45)

    # Adding title and labels
    plt.title(f'Weights vs Counts after sampling {n} times')
    plt.xlabel('Weights')
    plt.ylabel('Counts')

    # Ensuring layout is tight to accommodate rotated x-axis labels
    plt.tight_layout()

    plt.savefig(f'{outPrefix}.png', format='png')
    plt.savefig(f'{outPrefix}.pdf', format='pdf')

if __name__ == "__main__":
    main()