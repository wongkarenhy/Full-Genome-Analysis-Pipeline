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
from .tenxLargeSvDeletions import read10xlargeSVs




def tenxlargesvunknown(args):

    #load sample data
    sample_frame = read10xlargeSVs(args.samplepath, '<UNK>', True)

    #load parent data
    father_frame = read10xlargeSVs(args.fpath, '<UNK>', False)
    mother_frame = read10xlargeSVs(args.mpath, '<UNK>', False)

    #load reference data
    ref_frame = read10xlargeSVs(args.referencepath, '<UNK>', False)

    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()
    sample_start, sample_end, mother_start, mother_end, father_start, father_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), mother_frame.copy(),  mother_frame.copy(), father_frame.copy(), father_frame.copy(), ref_frame.copy(), ref_frame.copy()
    
    for df in [sample_start, father_start, mother_start, ref_start]: #create an interval for the start point
        df['Start'], df['End'], df['Chromosome'] = df.POS - 10000, df.POS + 10000, df['CHROM']

    for df in [sample_end, father_end, mother_end, ref_end]: #create an interval for the end point
        df['Start'], df['End'], df['Chromosome'] = df.END - 10000, df.END + 10000, df['CHROM']

    #overlap start and end points with genes separately  
    sample_frame =  geneOverlapINVBND(args, sample_start, sample_end, sample_frame)

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlapINVBND(sample_start, sample_end, ref_start, ref_end, sample_frame)
    if filtered_sample_frame.empty:
        filtered_sample_frame.to_csv(args.outputdirectory + '/confident_set/' + args.sampleID + '_10x_unknown_largeSV.txt', sep='\t', index=False)
        return None

    #add column based on overlap with parent
    if args.singleton:
        calls = checkParentsOverlapINVBND(sample_start, father_start, mother_start, sample_end, father_end, mother_end, filtered_sample_frame)
        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL', 'FILTER_PASS', 'END', 'SVLEN', 'Gene', 'OMIM_syndrome', 'Gene2', 'OMIM_syndrome2', 'Found_in_Father', 'Found_in_Mother']
    else:
        calls = filtered_sample_frame.rename(columns={'Name':'Gene', 'Name2':'Gene2', 'Score': 'OMIM_syndrome', 'Score2': 'OMIM_syndrome2'})
        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL', 'FILTER_PASS', 'END', 'SVLEN', 'Gene', 'OMIM_syndrome', 'Gene2', 'OMIM_syndrome2']

    # Write final output

    calls = calls[cols]
    calls.to_csv(args.outputdirectory + '/confident_set/' + args.sampleID + '_10x_unknown_largeSV.txt', sep='\t', index = False)




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

    tenxlargesvunknown(args)




if __name__=="__main__":
    main()
