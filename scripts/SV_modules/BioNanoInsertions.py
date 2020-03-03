#!/usr/bin/env python3.6

from .BioNanoDeletions import readsmap, reciprocal_overlap, checkRefOverlap, checkParentsOverlap, exonOverlap, overlap_length
import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse



def BN_insertion(args):

    #loadsample
    sample_frame = readsmap(args.samplepath, args, 'insertion')

    #loadparent
    mother_frame = readsmap(args.mpath, args, 'insertion')
    father_frame = readsmap(args.fpath, args, 'insertion')

    #load reference
    ref_frame = readsmap(args.referencepath, args, 'insertion')

    # Actual fun
    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']
        
    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)

    #add column based on overlap with parents
    if not args.singleton:
        df = checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame)
    else:
        df = filtered_sample_frame

    #describe exon overlap
    df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']
    exon_calls = exonOverlap(args, df)

    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_Bionano_insertions.txt', sep='\t', index = False)


def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-c", "--confidence", help="Give the confidence level cutoff for the sample here", dest="confidence", type=str, default=0.5)
    parser.add_argument("-e", "--exons", help="Give the file with exons intervals, names, and phenotypes here", dest="exons", type=str, required=True)
    parser.add_argument("-g", "--genelist", help="Primary genelist with scores", dest="genelist", type=str)
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()

    # Actual function
    BN_insertion(args)

if __name__=="__main__":
    main()
