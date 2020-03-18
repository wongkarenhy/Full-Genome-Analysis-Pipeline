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
    df = df[df['ALT_1'].notna()]
    df = df.loc[df['ALT_1'].str.contains("DUP")]
    df.reset_index(inplace=True, drop=True)
    df['CHROM'] = df['CHROM'].map(lambda x: x.lstrip('chr'))

    return(df)


def tenxlargesvduplications(args):

    #load sample data
    sample_frame = read10xlargeSVs(args.samplepath, 'DUP')

    if not args.singleton:
        #load parent data
        father_frame = read10xlargeSVs(args.fpath, 'DUP')
        mother_frame = read10xlargeSVs(args.mpath, 'DUP')
    else:
        mother_frame = pd.DataFrame(columns=['POS', 'END', 'CHROM'])
        father_frame = pd.DataFrame(columns=['POS', 'END', 'CHROM'])


    #load reference data
    ref_frame = read10xlargeSVs(args.referencepath, 'DUP')



    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)

    #add column based on overlap with parents
    if not args.singleton:
        df = checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame)
    else:
        df = filtered_sample_frame

    #describe exon overlap
    df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']
    exon_calls = exonOverlap(args, df)

    # Write final output
    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_10x_duplications_largeSV_exons.txt', sep='\t', index = False)

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
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()

    # Actual function
    tenxlargesvduplications(args)


if __name__=="__main__":
    main()

