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


def read10xlargeSVs(input, sv_type, qual_filter):

    df = allel.vcf_to_dataframe(input, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END', 'variants/SVLEN'])
    if qual_filter:
        scores_cutoff = np.mean(df.QUAL) + 1*np.std(df.QUAL)
        df = df.loc[df['QUAL'] > scores_cutoff]
    df = df.loc[df['ALT_1']==sv_type]
    df.reset_index(inplace=True, drop=True)
    df['CHROM'] = df['CHROM'].map(lambda x: x.lstrip('chr'))

    return(df)



def tenxlargesvdeletions(args):

    # Load sample data
    sample_frame = read10xlargeSVs(args.samplepath, '<DEL>', False)

    # Load parent data
    if args.type == 'trio' or args.father_duo:
        father_frame = read10xlargeSVs(args.fpath, '<DEL>', False)
    else:
        father_frame = pd.DataFrame(columns=['POS', 'END', 'CHROM'])

    if args.type == 'trio' or args.mother_duo:
        mother_frame = read10xlargeSVs(args.mpath, '<DEL>', False)
    else:
        mother_frame = pd.DataFrame(columns=['POS', 'END', 'CHROM'])

    # Load reference data
    ref_frame = read10xlargeSVs(args.referencepath, '<DEL>', False)

    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']

    # Remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)

    # Add column based on overlap with parents
    filtered_sample_frame = checkParentsOverlap(sample_copy, father_copy, filtered_sample_frame, args, 'Found_in_Father')
    filtered_sample_frame = checkParentsOverlap(sample_copy, mother_copy, filtered_sample_frame, args, 'Found_in_Mother')

    # Describe exon overlap
    filtered_sample_frame['Start'], filtered_sample_frame['End'], filtered_sample_frame['Chromosome'] = filtered_sample_frame.POS, filtered_sample_frame.END, filtered_sample_frame['CHROM']
    exon_calls = exonOverlap(args, filtered_sample_frame)

    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_10x_deletions_largeSV_exons.txt', sep='\t', index = False)

    return filtered_sample_frame, exon_calls

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
    parser.add_argument("-t", "--type", help="Specify whether this is a trio, duo, or singleton case", dest="type", type=str)
    parser.add_argument("-F", help="Set this flag if this is a duo case AND only father is sequenced", dest="father_duo", action='store_true')
    parser.add_argument("-M", help="Set this flag if this is a duo case AND only mother is sequenced", dest="mother_duo", action='store_true')
    args = parser.parse_args()

    # Actual function
    tenxlargesvdeletions(args)


if __name__=="__main__":
    main()

