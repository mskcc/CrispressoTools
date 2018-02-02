#!/usr/bin/env python2.7

from Bio import SeqIO
from collections import Counter
import gzip
import argparse
import sys
import os
import numpy
import logging
logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info

__version__ = "0.9.4"

nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})

def reverse_complement(seq):
        return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference(set(['A','T','C','G','N'])))

def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise Exception('I cannot open the file: '+filename)

def readFASTQ(fastq_filename):
    if fastq_filename.endswith('.gz'):
        fastq_handle=gzip.open(fastq_filename)
    else:
        fastq_handle=open(fastq_filename)

    for record in SeqIO.parse(fastq_handle, "fastq"):
        yield record

def find_tag_on_reference(tag_len, read1, read2, ref_sequence, reverse_complement):
    print "Processing forward reads"
    foundR1tag=False
    foundR2tag=False

    results = []

    end_R1tag = tag_len
    end_R2tag = tag_len
    
    start_R1tag=0
    start_R2tag=0

    tagR1_sequence = read1[-end_R1tag:]
    tagR2_sequence = read2[-end_R2tag:]


    pos = ref_sequence.find(str(tagR1_sequence.seq))
    pos2 = reverse_complement.find(str(tagR2_sequence.seq))

    if pos>-1:
        foundR1tag=True
        print "found R1tag"
        results.append(pos)
        results.append(tagR1_sequence)
        results.append(end_R1tag)

    if pos2>-1:
        foundR2tag=True
        print "found R2tag"
        results.append(pos2)
        results.append(tagR2_sequence)
        results.append(end_R2tag)

    if not foundR1tag and not pos>-1:
        i = 1
        max_iterations = len(ref_sequence)%tag_len
        
        while i<=max_iterations and not foundR1tag:
            print "finding R1tag"
            i +=1
            end_R1tag+=tag_len
            start_R1tag+=tag_len
            tagR1_sequence=read1[-end_R1tag: -start_R1tag]
            pos = ref_sequence.find(str(tagR1_sequence.seq))

            if pos>-1:
                print "Found with HOPPING"
                foundR1tag = True
                results.append(pos)
                results.append(tagR1_sequence)
                results.append(end_R1tag)

    if not foundR2tag and not pos2>-1:
        j = 1
        max_iterations = len(ref_sequence)%tag_len

        while j<=max_iterations and not foundR2tag:
            print "finding R2tag" 
            j +=1
            end_R2tag+=tag_len
            start_R2tag+=tag_len
            tagR2_sequence=read2[-end_R2tag: -start_R2tag]
            pos2 = reverse_complement.find(str(tagR2_sequence.seq))

            if pos2>-1:
                print "Found with HOPPING"
                foundR2tag = True

                results.append(pos2)
                results.append(tagR2_sequence)
                results.append(end_R2tag)


    if not pos>-1 or not pos2>-1:
        results.append(pos)
        results.append(tagR1_sequence)
        results.append(end_R1tag)
        results.append(pos2)
        results.append(tagR2_sequence)
        results.append(end_R2tag)

   # print results
    return results 

def find_tag_on_revCompliment(tag_len, read1, read2, ref_sequence, reverse_complement):
    print "Processing reverse reads"
    foundR1tag=False
    foundR2tag=False

    results = []

    end_R1tag = tag_len
    end_R2tag = tag_len
    
    start_R1tag=0
    start_R2tag=0

    tagR1_sequence = read1[-end_R1tag:]
    tagR2_sequence = read2[-end_R2tag:]


    posReverseComp = reverse_complement.find(str(tagR1_sequence.seq))
    pos2ReverseComp = ref_sequence.find(str(tagR2_sequence.seq))

    if posReverseComp>-1:
        print "found rev R1tag"
        foundR1tag=True
        results.append(posReverseComp)
        results.append(tagR1_sequence)
        results.append(end_R1tag)

    if pos2ReverseComp>-1:
        foundR2tag=True
        print "found rev R2tag"
        results.append(pos2ReverseComp)
        results.append(tagR2_sequence)
        results.append(end_R2tag)

    if not foundR1tag and not posReverseComp>-1:
        i = 1
        max_iterations = len(ref_sequence)%tag_len

        while i<=max_iterations and not foundR1tag:
            print "finding rev  R1tag"
            i +=1
            end_R1tag+=tag_len
            start_R1tag+=tag_len
            tagR1_sequence=read1[-end_R1tag: -start_R1tag]
            posReverseComp = reverse_complement.find(str(tagR1_sequence.seq))

            if posReverseComp>-1:
                print "Found with HOPPING"
                foundR1tag = True
                results.append(posReverseComp)
                results.append(tagR1_sequence)
                results.append(end_R1tag)

    if not foundR2tag and not pos2ReverseComp>-1:
        j = 1
        max_iterations = len(ref_sequence)%tag_len
        while j<=max_iterations and not foundR2tag:
            print "finding rev  R2tag"
            j +=1
            end_R2tag+=tag_len
            start_R2tag+=tag_len
            tagR2_sequence=read2[-end_R2tag: -start_R2tag]
            pos2ReverseComp = ref_sequence.find(str(tagR2_sequence.seq))

            if pos2ReverseComp>-1:
                print "Found with HOPPING"
                foundR2tag = True
                results.append(pos2ReverseComp)
                results.append(tagR2_sequence)
                results.append(end_R2tag)

    if not posReverseComp>-1 or not pos2ReverseComp>-1:

        results.append(posReverseComp)
        results.append(tagR1_sequence)
        results.append(end_R1tag)
        results.append(pos2ReverseComp)
        results.append(tagR2_sequence)
        results.append(end_R2tag)
        
    #print results
    return results

def main():
    # try:
        print 'Version %s\n' % __version__

        parser = argparse.ArgumentParser(
            description='fixNonOverlappingReads Parameters',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-r1','--fastq_r1', type=str,
            help='First fastq file', required=True,default='Fastq filename' )
        parser.add_argument('-r2','--fastq_r2', type=str,
            help='Second fastq file for paired end reads', required=True,default='')
        parser.add_argument('-a','--amplicon_seq', type=str,
            help='Amplicon Sequence', required=True)
        parser.add_argument('-t','--tag_len', type=int,
            help='Tag (anchor) length', default=10)
        args = parser.parse_args()

        #check files
        check_file(args.fastq_r1)
        if args.fastq_r2:
             check_file(args.fastq_r2)

        #amplicon sequence check
        #make evetything uppercase!
        args.amplicon_seq=args.amplicon_seq.strip().upper()
        revComp_amplicon=reverse_complement(args.amplicon_seq)
        wrong_nt=find_wrong_nt(args.amplicon_seq)
        if wrong_nt:
            raise NTException('The amplicon sequence contains wrong characters:%s' % ' '.join(wrong_nt))

        len_amplicon=len(args.amplicon_seq)

        print "Amplicon length =",len_amplicon

        NUM_SEQS_TO_SCAN=10000

        newR1File=os.path.basename(args.fastq_r1).replace(".fastq.gz","__CP_PAD.fastq.gz")
        fpR1=gzip.open(newR1File,"w")
        newR2File=os.path.basename(args.fastq_r2).replace(".fastq.gz","__CP_PAD.fastq.gz")
        fpR2=gzip.open(newR2File,"w")
        readsNotPaddedR1=os.path.basename(args.fastq_r1).replace(".fastq.gz","__CP_NotPadded.fastq.gz")
        notpaddedR1= gzip.open(readsNotPaddedR1, "w")
        readsNotPaddedR2=os.path.basename(args.fastq_r2).replace(".fastq.gz","__CP_NotPadded.fastq.gz")
        notpaddedR2=gzip.open(readsNotPaddedR2, "w")

        print newR1File, newR2File, readsNotPaddedR1, readsNotPaddedR2
        total_reads_processed=0
        total_reads_padded=0
        notFoundCount = 0
        
        for i,rr in enumerate(zip(readFASTQ(args.fastq_r1),readFASTQ(args.fastq_r2))):
            r1,r2=rr
            

            total_reads_processed +=2

            tag_forward_reads = find_tag_on_reference(args.tag_len, r1, r2, args.amplicon_seq, revComp_amplicon)
            #tag_reverse_reads = find_tag_on_revCompliment(args.tag_len, r1, r2, args.amplicon_seq, revComp_amplicon)
            # tagR1=r1[-args.tag_len:]
            # tagR2=r2[-args.tag_len:]

            tagR1_forward = tag_forward_reads[1]
            tagR2_forward = tag_forward_reads[4]

            tagR1_forward_distance = tag_forward_reads[2]
            tagR2_forward_distance = tag_forward_reads[5]

            pos=tag_forward_reads[0]
            pos2=tag_forward_reads[3]
          
	    posReverseStrand=-1
            pos2ReverseStrand=-1
	    if not pos>-1 and not pos2>-1:
                tag_reverse_reads = find_tag_on_revCompliment(args.tag_len, r1, r2, args.amplicon_seq, revComp_amplicon) 
                tagR1_rev = tag_reverse_reads[1]
                tagR2_rev = tag_reverse_reads[4]

                tagR1_rev_distance = tag_reverse_reads[2]
                tagR2_rev_distance = tag_reverse_reads[5]

                posReverseStrand=tag_reverse_reads[0]
                pos2ReverseStrand=tag_reverse_reads[3]

            padLen=150-len(r1)
            # pos=args.amplicon_seq.find(str(tagR1.seq))
            # pos2=revComp_amplicon.find(str(tagR2.seq))
            # posReverseStrand=revComp_amplicon.find(str(tagR1.seq))
            # pos2ReverseStrand=args.amplicon_seq.find(str(tagR2.seq))
            
            print pos, pos2, posReverseStrand, pos2ReverseStrand

            if pos>-1 and pos2>-1:
                #print pos, pos2

                qVals=[x for x in tagR1_forward.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r1.letter_annotations["phred_quality"]])
                extendedRead1= str(r1.seq)+args.amplicon_seq[(pos+tagR1_forward_distance):(pos+padLen+tagR1_forward_distance)]

                print >>fpR1, r1.description
                print >>fpR1, extendedRead1
                print >>fpR1, "+"
                print >>fpR1, (origQ+extendQ)[:len(extendedRead1)]

                qVals=[x for x in tagR2_forward.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r2.letter_annotations["phred_quality"]])
                extendedRead2= str(r2.seq)+revComp_amplicon[(pos2+tagR2_forward_distance):(pos2+padLen+tagR2_forward_distance)]

                print >>fpR2, r2.description
                print >>fpR2, extendedRead2
                print >>fpR2, "+"
                print >>fpR2, (origQ+extendQ)[:len(extendedRead2)]

                total_reads_padded+=2

            elif  posReverseStrand>-1 and pos2ReverseStrand>-1:
                #print str(posReverseStrand) + " Reverse" , str(pos2ReverseStrand) + " Reverse"

                qVals=[x for x in tagR1_rev.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r1.letter_annotations["phred_quality"]])
                extendedRead1=str(r1.seq)+revComp_amplicon[(posReverseStrand+tagR1_rev_distance):(posReverseStrand+padLen+tagR1_rev_distance)]

                print >>fpR1, r1.description
                print >>fpR1, extendedRead1
                print >>fpR1, "+"
                print >>fpR1, (origQ+extendQ)[:len(extendedRead1)]

                qVals=[x for x in tagR2_rev.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r2.letter_annotations["phred_quality"]])
                extendedRead2=str(r2.seq)+args.amplicon_seq[(pos2ReverseStrand+tagR1_rev_distance):(pos2ReverseStrand+padLen+tagR1_rev_distance)]

                print >>fpR2, r2.description
                print >>fpR2, extendedRead2
                print >>fpR2, "+"
                print >>fpR2, (origQ+extendQ)[:len(extendedRead2)]

                total_reads_padded+=2

           # if (not pos>-1 or not pos2>-1) and (not posReverseStrand>-1 or pos2ReverseStrand>-1):
            else:    
                print >>notpaddedR1, r1.description
                print >>notpaddedR1, r1.seq
                print >>notpaddedR1, "+"
                print >>notpaddedR1, "".join([chr(x+33) for x in r1.letter_annotations["phred_quality"]])
                
                print >>notpaddedR2, r2.description
                print >>notpaddedR2, r2.seq
                print >>notpaddedR2, "+"
                print >>notpaddedR2, "".join([chr(x+33) for x in r2.letter_annotations["phred_quality"]])

                notFoundCount+=2

            if i>NUM_SEQS_TO_SCAN:

                break
                 
        print "Total reads processed, " + str(total_reads_processed)
        print "Total reads padded, " + str(total_reads_padded)
        print "Total reads eliminated in this process, " + str(notFoundCount)

        fpR1.close()
        fpR2.close()

    # except Exception as e:
    #      error('Unexpected error, please check your input.\n\nERROR: %s' % e)
    #      sys.exit(-1)

if __name__=="__main__":
    main()

