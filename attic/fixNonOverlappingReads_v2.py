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

__version__ = "0.9.1"

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

        NUM_SEQS_TO_SCAN=100

        newR1File=os.path.basename(args.fastq_r1).replace(".fastq.gz","__CP_PAD.fastq.gz")
        fpR1=gzip.open(newR1File,"w")
        newR2File=os.path.basename(args.fastq_r2).replace(".fastq.gz","__CP_PAD.fastq.gz")
        fpR2=gzip.open(newR2File,"w")

        print newR1File, newR2File

        for i,rr in enumerate(zip(readFASTQ(args.fastq_r1),readFASTQ(args.fastq_r2))):
            r1,r2=rr
            tagR1=r1[-args.tag_len:]
            tagR2=r2[-args.tag_len:]
            padLen=150-len(r1)
            pos=args.amplicon_seq.find(str(tagR1.seq))
            pos2=revComp_amplicon.find(str(tagR2.seq))
            posReverseStrand=revComp_amplicon.find(str(tagR1.seq))
            pos2ReverseStrand=args.amplicon_seq.find(str(tagR2.seq)) 
            if pos>-1 and pos2>-1:
                print pos, pos2
                qVals=[x for x in tagR1.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r1.letter_annotations["phred_quality"]])
                extendedRead1= str(r1.seq)+args.amplicon_seq[(pos+args.tag_len):(pos+padLen+args.tag_len)]

                print >>fpR1, r1.description
                print >>fpR1, extendedRead1
                print >>fpR1, "+"
                print >>fpR1, (origQ+extendQ)[:len(extendedRead1)]

                qVals=[x for x in tagR2.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r2.letter_annotations["phred_quality"]])
                extendedRead2= str(r2.seq)+revComp_amplicon[(pos2+args.tag_len):(pos2+padLen+args.tag_len)]

                print >>fpR2, r2.description
                print >>fpR2, extendedRead2
                print >>fpR2, "+"
                print >>fpR2, (origQ+extendQ)[:len(extendedRead2)]


            if  pos==-1 and pos2==-1 and posReverseStrand>-1 and pos2ReverseStrand>-1:
                print str(posReverseStrand) + " Reverse" , str(pos2ReverseStrand) + " Reverse"

                qVals=[x for x in tagR1.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r1.letter_annotations["phred_quality"]])
                extendedRead1=str(r1.seq)+revComp_amplicon[(posReverseStrand+args.tag_len):(posReverseStrand+padLen+args.tag_len)]

                print >>fpR1, r1.description
                print >>fpR1, extendedRead1
                print >>fpR1, "+"
                print >>fpR1, (origQ+extendQ)[:len(extendedRead1)]

                qVals=[x for x in tagR2.letter_annotations["phred_quality"]]
                meanQ=int(numpy.mean(qVals))
                extendQ=chr(meanQ+33)*padLen
                origQ="".join([chr(x+33) for x in r2.letter_annotations["phred_quality"]])
                extendedRead2=str(r2.seq)+args.amplicon_seq[(pos2ReverseStrand+args.tag_len):(pos2ReverseStrand+padLen+args.tag_len)]

                print >>fpR2, r2.description
                print >>fpR2, extendedRead2
                print >>fpR2, "+"
                print >>fpR2, (origQ+extendQ)[:len(extendedRead2)]


            # if i>NUM_SEQS_TO_SCAN:
            #     break

        fpR1.close()
        fpR2.close()

    # except Exception as e:
    #      error('Unexpected error, please check your input.\n\nERROR: %s' % e)
    #      sys.exit(-1)

if __name__=="__main__":
    main()

