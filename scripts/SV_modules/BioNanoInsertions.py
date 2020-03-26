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

    # Load Sample
    sample_frame = readsmap(args.samplepath, args, 'insertion')

    # Load parent data
    if args.type == 'trio' or args.father_duo:
        father_frame = readsmap(args.fpath, args, 'insertion')
    else:
        father_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1'])

    if args.type == 'trio' or args.mother_duo:
        mother_frame = readsmap(args.mpath, args, 'insertion')
    else:
        mother_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1'])

    # Load reference
    ref_frame = readsmap(args.referencepath, args, 'insertion')

    # Actual fun
    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']
        
    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)

    #add column based on overlap with parents
    filtered_sample_frame = checkParentsOverlap(sample_copy, father_copy, filtered_sample_frame, args, 'Found_in_Father')
    filtered_sample_frame = checkParentsOverlap(sample_copy, mother_copy, filtered_sample_frame, args, 'Found_in_Mother')

    #describe exon overlap
    filtered_sample_frame['Start'], filtered_sample_frame['End'], filtered_sample_frame['Chromosome'] = filtered_sample_frame.RefStartPos, filtered_sample_frame.RefEndPos, filtered_sample_frame['RefcontigID1']
    exon_calls = exonOverlap(args, filtered_sample_frame)

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
    parser.add_argument("-t", "--type", help="Specify whether this is a trio, duo, or singleton case", dest="type", type=str)
    parser.add_argument("-F", help="Set this flag if this is a duo case AND only father is sequenced", dest="father_duo", action='store_true')
    parser.add_argument("-M", help="Set this flag if this is a duo case AND only mother is sequenced", dest="mother_duo", action='store_true')
    args = parser.parse_args()

    # Actual function
    BN_insertion(args)

if __name__=="__main__":
    main()
