CC = g++
COFLAGS = -Wall -O3 -c -I.
PROGRAMS = csem csem-bam2wig extractFromEland csem-bam-processor

all : $(PROGRAMS)

sam/libbam.a :
	cd sam ; ${MAKE} all

BamAlignment.h : sam/bam.h sam/sam.h utils.h my_assert.h

SamParser.h : sam/bam.h sam/sam.h utils.h my_assert.h BamAlignment.h

ChrMap.h : sam/bam.h utils.h

sam_csem_aux.h : sam/bam.h

BamWriter.h : sam/bam.h sam/sam.h utils.h my_assert.h sam_csem_aux.h BamAlignment.h

Alignment.h : utils.h

ArrayScan.h : utils.h

Chromosome.h : utils.h Alignment.h ArrayScan.h

ChromTable.h : utils.h my_assert.h ChrMap.h Alignment.h Chromosome.h PThreadWrapper.h

csem.o : sam/bam.h sam/sam.h utils.h my_assert.h BamAlignment.h SamParser.h ChrMap.h BamWriter.h Alignment.h ArrayScan.h Chromosome.h ChromTable.h PThreadWrapper.h csem.cpp
	$(CC) $(COFLAGS) -ffast-math csem.cpp 

csem : csem.o sam/libbam.a
	$(CC) -o $@ csem.o sam/libbam.a -lz -lpthread

wiggle.cpp : utils.h wiggle.h

wiggle.o : sam/bam.h sam/sam.h utils.h wiggle.h wiggle.cpp
	$(CC) $(COFLAGS) wiggle.cpp

bam2wig.o : wiggle.h bam2wig.cpp
	$(CC) $(COFLAGS) bam2wig.cpp

csem-bam2wig : wiggle.o bam2wig.o sam/libbam.a
	$(CC) -o $@ wiggle.o bam2wig.o sam/libbam.a -lz

extractFromEland.o : extractFromEland.cpp
	$(CC) $(COFLAGS) extractFromEland.cpp

extractFromEland : extractFromEland.o
	$(CC) -o $@ extractFromEland.o

sampling.h : boost/random.hpp

bamProcessor.o : sam/bam.h sam/sam.h utils.h my_assert.h sampling.h BamAlignment.h SamParser.h BamWriter.h bamProcessor.cpp
	$(CC) $(COFLAGS) bamProcessor.cpp

csem-bam-processor : bamProcessor.o sam/libbam.a
	$(CC) -o $@ bamProcessor.o sam/libbam.a -lz 

clean :
	rm -f *.o *~ $(PROGRAMS)
	cd sam ; ${MAKE} clean
