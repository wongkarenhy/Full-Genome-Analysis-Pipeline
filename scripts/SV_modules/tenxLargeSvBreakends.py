#!/usr/bin/env python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import allel
import argparse
from .tenxLargeSvInversions import geneOverlapINVBND, checkRefOverlapINVBND, checkParentsOverlapINVBND



def read_breakend(path, qual_filter):
    sample_frame = allel.vcf_to_dataframe(path, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF',
                                                        'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS',
                                                        'variants/MATEID', 'variants/SVTYPE'])
    if qual_filter:
        scores_cutoff = np.mean(sample_frame.QUAL) + 1*np.std(sample_frame.QUAL)
        sample_frame = sample_frame.loc[sample_frame['QUAL'] > scores_cutoff]
    df = sample_frame.loc[sample_frame['SVTYPE'].str.contains("BND", na=False)]
    df = pd.merge(df, df[['ID', 'MATEID', 'CHROM', 'POS']], left_on="ID", right_on="MATEID")
    df = df[~df.MATEID_x.str.endswith("_1")]
    df = df[~df.MATEID_x.str.endswith("_3")]
    df = df.rename(columns={'ID_x': 'ID'})
    df['CHROM_x'] = df['CHROM_x'].map(lambda x: x.lstrip('chr'))
    df['CHROM_y'] = df['CHROM_y'].map(lambda x: x.lstrip('chr'))

    return df

def tenxlargesvbreakends(args):

    sample_frame = read_breakend(args.samplepath, True)
    mother_frame = read_breakend(args.mpath, False)
    father_frame = read_breakend(args.fpath, False)
    ref_frame = read_breakend(args.referencepath, False)

    sample_start, sample_end, mother_start, mother_end, father_start, father_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), mother_frame.copy(),  mother_frame.copy(), father_frame.copy(), father_frame.copy(), ref_frame.copy(), ref_frame.copy()
    for df in [sample_start, father_start, mother_start, ref_start]: #create an interval for the breakend start point
        df['Start'], df['End'], df['Chromosome'] = df.POS_x - 10000, df.POS_x + 10000, df['CHROM_x']

    for df in [sample_end, father_end, mother_end, ref_end]: #create an interval for the breakend end point
        df['Start'], df['End'], df['Chromosome'] = df.POS_y - 10000, df.POS_y + 10000, df['CHROM_y']

    #overlap start and end points with genes separately
    sample_frame =  geneOverlapINVBND(args, sample_start, sample_end, sample_frame)

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlapINVBND(sample_start, sample_end, ref_start, ref_end, sample_frame)
    if filtered_sample_frame.empty:
        filtered_sample_frame.to_csv(args.outputdirectory + '/' + args.sampleID + '_10x_breakends_largeSV.txt', sep='\t', index=False)
        return None

    #add column based on overlap with parent
    if not args.singleton:
        calls = checkParentsOverlapINVBND(sample_start, father_start, mother_start, sample_end, father_end, mother_end, filtered_sample_frame)
        cols = ['CHROM_x', 'POS_x', 'ID', 'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL', 'FILTER_PASS', 'SVTYPE', 'MATEID_x', 'Gene', 'OMIM_syndrome', 'Gene2', 'OMIM_syndrome2', 'Found_in_Father', 'Found_in_Mother']
    else:
        calls = filtered_sample_frame.rename(columns={'Name':'Gene', 'Name2':'Gene2', 'Score': 'OMIM_syndrome', 'Score2': 'OMIM_syndrome2'})
    cols = ['CHROM_x', 'POS_x', 'ID', 'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL', 'FILTER_PASS', 'SVTYPE', 'MATEID_x','Gene', 'OMIM_syndrome', 'Gene2', 'OMIM_syndrome2']

    # Write final output
    calls = calls[cols]
    calls.to_csv(args.outputdirectory + '/confident_set/' + args.sampleID + '_10x_breakends_largeSV.txt', sep='\t', index = False)

    return None




def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-e", "--genes", help="Give the file with genes intervals, names, and phenotypes here", dest="genes", type=str, required=True)
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()

    tenxlargesvbreakends(args)

if __name__=="__main__":
    main()
