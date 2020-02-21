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

def tenxlargesvinversions(args):

    #load sample data
    sample_frame = read10xlargeSVs(args.samplepath, '<INV>')

    #load parent data
    father_frame = read10xlargeSVs(args.fpath, '<INV>')
    mother_frame = read10xlargeSVs(args.mpath, '<INV>')

    #load reference data
    ref_frame = read10xlargeSVs(args.referencepath, '<INV>')

    sample_start, sample_end, mother_start, mother_end, father_start, father_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), mother_frame.copy(), mother_frame.copy(), father_frame.copy(), father_frame.copy(), ref_frame.copy(), ref_frame.copy()

    for df in [sample_start, father_start, mother_start, ref_start]: #create an interval for the inversion start point
        df['Start'], df['End'], df['Chromosome'] = df.POS - 100000, df.POS + 100000, df['CHROM']
                        
    for df in [sample_end, father_end, mother_end, ref_end]: #create an interval for the inversion end point
        df['Start'], df['End'], df['Chromosome'] = df.END - 100000, df.END + 100000, df['CHROM']
    
    #overlap start and end points with exons separately  
    exon_frame = pr.read_bed(args.exons)              
    exon_start = PyRanges(sample_start).join(exon_frame[["Name", "Score"]]).drop(like="_b")
    exon_end = PyRanges(sample_end).join(exon_frame[["Name", "Score"]]).drop(like="_b")  

    if exon_start.df.empty and exon_end.df.empty:
        sample_frame['Name'] = sample_frame['Name2'] = sample_frame['Score'] = sample_frame['Score2'] ='None' 
    elif exon_start.df.empty:
        sample_frame = sample_frame.merge(exon_end.df.rename(columns={'Name': 'Name2', 'Score': 'Score2'}).filter(items=['ID', 'Name']).drop_duplicates(), on='ID', how="left")
        #sample_frame = exon_end.df.rename(columns = {'Name':'Name2', 'Score':'Score2'}).filter(items=['ID', 'Name']).drop_duplicates().merge(sample_frame, on=['ID'], how='right')
        sample_frame['Name'] = sample_frame['Score'] = 'None'
    elif exon_end.df.empty:
        sample_frame = sample_frame.merge(exon_start.df.filter(items=['ID', 'Name', 'Score']).drop_duplicates(), on='ID', how='left')
        #sample_frame = exon_start.df.filter(items=['ID', 'Name', 'Score']).drop_duplicates().merge(sample_frame, on=['ID'], how='right')
        sample_frame['Name2'] = sample_frame['Score2'] = 'None'
    else: 
        sample_frame = sample_frame.merge(exon_start.df.filter(items=['ID', 'Name', 'Score']).drop_duplicates(), on='ID', how='left')
        #sample_frame = exon_start.df.filter(items=['ID', 'Name', 'Score']).drop_duplicates().merge(sample_frame, on=['ID'], how='right')
        sample_frame = sample_frame.merge(exon_end.df.rename(columns={'Name': 'Name2', 'Score': 'Score2'}).filter(items=['ID', 'Name2', 'Score2']),on=['ID'], how='left')

    #remove anything that overlaps with the reference
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
    
    #add column based on overlap with parent
    denovo_start_f, denovo_start_m = PyRanges(sample_start).overlap(PyRanges(father_start)), PyRanges(sample_start).overlap(PyRanges(mother_start))
    denovo_end_f, denovo_end_m = PyRanges(sample_end).overlap(PyRanges(father_end)), PyRanges(sample_end).overlap(PyRanges(mother_end))
    denovo_f_frame = pd.concat([denovo_start_f.df, denovo_end_f.df], axis=0)
    denovo_m_frame = pd.concat([denovo_start_m.df, denovo_end_m.df], axis=0)


    if denovo_f_frame.empty:
        filtered_sample_frame["Found_in_Father"] = "False"
        filtered_sample_frame = filtered_sample_frame.rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()
    else:
        filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_f_frame, on=None, how='left', indicator='Found_in_Father')
        filtered_sample_frame['Found_in_Father'] = np.where(filtered_sample_frame.Found_in_Father == 'both', True, False)
        filtered_sample_frame['CHROM'] = filtered_sample_frame['CHROM'].map(lambda x: x.lstrip('chr'))
        filtered_sample_frame = filtered_sample_frame.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()

    if denovo_m_frame.empty:
        filtered_sample_frame["Found_in_Mother"] = "False"
        calls = filtered_sample_frame.rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()
    else:
        filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_m_frame, on=None, how='left', indicator='Found_in_Mother')
        filtered_sample_frame['Found_in_Mother'] = np.where(filtered_sample_frame.Found_in_Mother == 'both', True, False)
        filtered_sample_frame['CHROM'] = filtered_sample_frame['CHROM'].map(lambda x: x.lstrip('chr'))
        calls = filtered_sample_frame.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()

    # Write final output
    cols = ['Name', 'Name2', 'Phenotype', 'Phenotype2', 'Found_in_Father', 'CHROM' ,'POS' ,'ID' ,'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL', 'FILTER_PASS', 'END', 'Found_in_Mother']
    calls = calls[cols]
    calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_10x_inversions_largeSV_exons.txt', sep='\t', index = False)


def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-e", "--exons", help="Give the file with exons intervals, names, and phenotypes here", dest="exons", type=str, required=True)
    args = parser.parse_args()

    # Actual function
    tenxlargesvinversions(args)


if __name__=="__main__":
    main()

