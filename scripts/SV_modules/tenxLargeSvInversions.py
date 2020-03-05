#!/usr/bin/env python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import allel
import argparse
from .tenxLargeSvDeletions import read10xlargeSVs


def geneOverlapINVBND(args, sample_start,sample_end, sample_frame):

    gene_frame = pr.read_bed(args.genes)
    gene_start = PyRanges(sample_start).join(gene_frame[["Name", "Score"]]).drop(like="_b")
    gene_end = PyRanges(sample_end).join(gene_frame[["Name", "Score"]]).drop(like="_b")

    if gene_start.df.empty and gene_end.df.empty:
        sample_frame['Name'] = sample_frame['Name2'] = sample_frame['Score'] = sample_frame['Score2'] ='None'
    elif gene_start.df.empty:
        sample_frame = sample_frame.merge(gene_end.df.rename(columns={'Name': 'Name2', 'Score': 'Score2'}).filter(items=['ID', 'Name']).drop_duplicates(), on='ID', how="left")
        sample_frame['Name'] = sample_frame['Score'] = 'None'
    elif gene_end.df.empty:
        sample_frame = sample_frame.merge(gene_start.df.filter(items=['ID', 'Name', 'Score']).drop_duplicates(), on='ID', how='left')
        sample_frame['Name2'] = sample_frame['Score2'] = 'None'
    else:
        sample_frame = sample_frame.merge(gene_start.df.filter(items=['ID', 'Name', 'Score']).drop_duplicates(), on='ID', how='left')
        sample_frame = sample_frame.merge(gene_end.df.rename(columns={'Name': 'Name2', 'Score': 'Score2'}).filter(items=['ID', 'Name2', 'Score2']),on=['ID'], how='left')

    return sample_frame



def checkRefOverlapINVBND(sample_start, sample_end, ref_start, ref_end, sample_frame):

    overlap_start = PyRanges(sample_start).overlap(PyRanges(ref_start))
    overlap_end = PyRanges(sample_end).overlap(PyRanges(ref_end))
    if overlap_start.df.empty and overlap_end.df.empty:
        filtered_sample_frame = sample_frame
    else:
        overlap_frame = overlap_start.df.merge(overlap_end.df, on=['ID'])
        if overlap_frame.empty:
            filtered_sample_frame = sample_frame
        else:
            common = sample_frame.merge(overlap_frame,on=['ID'])
            filtered_sample_frame = sample_frame[(~sample_frame.ID.isin(common.ID))]

    return filtered_sample_frame



def checkParentsOverlapINVBND(sample_start, father_start, mother_start, sample_end, father_end, mother_end, filtered_sample_frame):

    denovo_start_f, denovo_start_m = PyRanges(sample_start).overlap(PyRanges(father_start)), PyRanges(sample_start).overlap(PyRanges(mother_start))
    denovo_end_f, denovo_end_m = PyRanges(sample_end).overlap(PyRanges(father_end)), PyRanges(sample_end).overlap(PyRanges(mother_end))
    #print(denovo_start_m)
    #print(denovo_end_m)
    denovo_f_frame = pd.merge(denovo_start_f.df, denovo_end_f.df['ID'], on=['ID']).drop_duplicates()
    denovo_m_frame = pd.merge(denovo_start_m.df, denovo_end_m.df['ID'], on=['ID']).drop_duplicates()

    if denovo_f_frame.empty:
        filtered_sample_frame["Found_in_Father"] = "False"
    if denovo_m_frame.empty:
        filtered_sample_frame["Found_in_Mother"] = "False"
    elif not (denovo_f_frame.empty and denovo_m_frame.empty):
        f_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_f_frame, on=None, how='left',
                                           indicator='Found_in_Father')
        f_filtered_sample_frame['Found_in_Father'] = np.where(f_filtered_sample_frame.Found_in_Father == 'both', True,False)
        f_filtered_sample_frame = f_filtered_sample_frame.drop_duplicates().reset_index(drop=True)
        m_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_m_frame, on=None, how='left',
                                           indicator='Found_in_Mother')
        m_filtered_sample_frame['Found_in_Mother'] = np.where(m_filtered_sample_frame.Found_in_Mother == 'both', True,False)
        m_filtered_sample_frame = m_filtered_sample_frame.drop_duplicates().reset_index(drop=True)
        f_filtered_sample_frame['Found_in_Mother'] = m_filtered_sample_frame['Found_in_Mother']
        filtered_sample_frame = f_filtered_sample_frame

    calls = filtered_sample_frame.drop(columns=['Chromosome', 'Start', 'End']).rename(columns={'Name':'Gene', 'Name2':'Gene2', 'Score': 'OMIM_syndrome', 'Score2': 'OMIM_syndrome2'}).drop_duplicates()

    return calls



def tenxlargesvinversions(args):

    #load sample data
    sample_frame = read10xlargeSVs(args.samplepath, '<INV>', False)

    #load parent data
    father_frame = read10xlargeSVs(args.fpath, '<INV>', False)
    mother_frame = read10xlargeSVs(args.mpath, '<INV>', False)

    #load reference data
    ref_frame = read10xlargeSVs(args.referencepath, '<INV>', False)

    sample_start, sample_end, mother_start, mother_end, father_start, father_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), mother_frame.copy(), mother_frame.copy(), father_frame.copy(), father_frame.copy(), ref_frame.copy(), ref_frame.copy()

    for df in [sample_start, father_start, mother_start, ref_start]: #create an interval for the inversion start point
        df['Start'], df['End'], df['Chromosome'] = df.POS - 10000, df.POS + 10000, df['CHROM']
                        
    for df in [sample_end, father_end, mother_end, ref_end]: #create an interval for the inversion end point
        df['Start'], df['End'], df['Chromosome'] = df.END - 10000, df.END + 10000, df['CHROM']

    # #overlap start and end points with genes separately
    sample_frame =  geneOverlapINVBND(args, sample_start, sample_end, sample_frame)

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlapINVBND(sample_start, sample_end, ref_start, ref_end, sample_frame)
    
    #add column based on overlap with parent
    if not args.singleton:
        calls = checkParentsOverlapINVBND(sample_start, father_start, mother_start, sample_end, father_end, mother_end, filtered_sample_frame)
        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL', 'FILTER_PASS', 'END', 'SVLEN', 'Gene', 'OMIM_syndrome', 'Gene2', 'OMIM_syndrome2', 'Found_in_Father', 'Found_in_Mother']

    else:
        calls = filtered_sample_frame.rename(columns={'Name':'Gene', 'Name2':'Gene2', 'Score': 'OMIM_syndrome', 'Score2': 'OMIM_syndrome2'})
        cols = ['CHROM' ,'POS' ,'ID' ,'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL', 'FILTER_PASS', 'END', 'SVLEN','Gene', 'OMIM_syndrome', 'Gene2', 'OMIM_syndrome2']


    # Write final output
    calls = calls[cols].drop_duplicates()
    calls.to_csv(args.outputdirectory + '/confident_set/' + args.sampleID + '_10x_inversions_largeSV.txt', sep='\t', index = False)


def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-g", "--genes", help="Give the file with gene-level intervals, names, and phenotypes here", dest="genes", type=str, required=True)
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()

    # Actual function
    tenxlargesvinversions(args)


if __name__=="__main__":
    main()

