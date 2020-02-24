#!/usr/bin/env python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import allel
import argparse
from .tenxDeletions import checkRefOverlap, checkParentsOverlap, exonOverlap, overlap_length, reciprocal_overlap


def read10xlargeSVs(input, sv_type):

    df = allel.vcf_to_dataframe(input, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END', 'variants/SVLEN'])
    # scores_cutoff = np.mean(df.QUAL) - 2*np.std(df.QUAL)
    # df = df.loc[df['FILTER_PASS']==True]
    # df = df.loc[df['QUAL']>scores_cutoff]
    df = df.loc[df['ALT_1']==sv_type]
    df.reset_index(inplace=True, drop=True)

    return(df)

def tenxlargesvdeletions(args):
                        

    #load sample data
    sample_frame = read10xlargeSVs(args.samplepath, '<DEL>')

    #load parent data
    father_frame = read10xlargeSVs(args.fpath, '<DEL>')
    mother_frame = read10xlargeSVs(args.mpath, '<DEL>')

    #load reference data
    ref_frame = read10xlargeSVs(args.referencepath, '<DEL>')


    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)

    #add column based on overlap with parents 
    df = checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame)
    #df.to_csv(args.outputdirectory + '/' + args.sampleID + '_10xLargeSVDeletions_cytobands.txt', sep='\t', index = False)

    #describe exon overlap
    exon_calls = exonOverlap(args, df)

    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_10x_deletions_largeSV_exons.txt', sep='\t', index = False)

    return df, exon_calls

def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-e", "--exons", help="Give the file with exons intervals, names, and phenotypes here", dest="exons", type=str, required=True)
    parser.add_argument("-g", "--genelist", help="Primary genelist with scores", dest="genelist", type=str)
    args = parser.parse_args()

    # Actual function
    tenxlargesvdeletions(args)


if __name__=="__main__":
    main()

