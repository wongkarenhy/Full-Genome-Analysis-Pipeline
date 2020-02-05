#!/usr/bin/env python3

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import allel
import argparse

def tenxlargesvinversions():
    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-e", "--exons", help="Give the file with exons intervals, names, and phenotypes here", dest="exons", type=str, required=True)
    parser.add_argument("-y", "--cytobands", help="Give the BED file with cytoband intervals", dest="cytobands", type=str, required=True)
    args = parser.parse_args()

    #load sample data 
    sample_frame = allel.vcf_to_dataframe(args.samplepath, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END'])
    scores_cutoff = np.mean(sample_frame.QUAL) - 2*np.std(sample_frame.QUAL)
    sample_frame = sample_frame.loc[sample_frame['FILTER_PASS']==True].loc[sample_frame['QUAL']>scores_cutoff].loc[sample_frame['ALT_1']=='<INV>']

    #load parent data 
    father_frame = allel.vcf_to_dataframe(args.fpath, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END'])
    scores_cutoff = np.mean(father_frame.QUAL) - 2*np.std(father_frame.QUAL)
    father_frame = father_frame.loc[father_frame['FILTER_PASS']==True].loc[father_frame['QUAL']>scores_cutoff].loc[father_frame['ALT_1']=='<INV>']
    mother_frame = allel.vcf_to_dataframe(args.mpath, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END'])
    scores_cutoff = np.mean(mother_frame.QUAL) - 2*np.std(mother_frame.QUAL)
    mother_frame = mother_frame.loc[mother_frame['FILTER_PASS']==True].loc[mother_frame['QUAL']>scores_cutoff].loc[mother_frame['ALT_1']=='<INV>']
    parent_frame = pd.concat([mother_frame, father_frame])

    #load reference data
    ref_frame = allel.vcf_to_dataframe(args.referencepath, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END'])
    scores_cutoff = np.mean(ref_frame.QUAL) - 2*np.std(ref_frame.QUAL)
    ref_frame = ref_frame.loc[ref_frame['FILTER_PASS']==True]
    ref_frame = ref_frame.loc[ref_frame['QUAL']>scores_cutoff]
    ref_frame = ref_frame.loc[ref_frame['ALT_1']=='<INV>']

    sample_start, sample_end, parent_start, parent_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), parent_frame.copy(), parent_frame.copy(), ref_frame.copy(), ref_frame.copy()

    for df in [sample_start, parent_start, ref_start]: #create an interval for the inversion start point
        df['Start'], df['End'], df['Chromosome'] = df.POS - 100000, df.POS + 100000, df['CHROM']
                        
    for df in [sample_end, parent_end, ref_end]: #create an interval for the inversion end point
        df['Start'], df['End'], df['Chromosome'] = df.END - 100000, df.END + 100000, df['CHROM']
    
    #overlap start and end points with exons separately  
    exon_frame = pr.read_bed(args.exons)              
    exon_start = PyRanges(sample_start).join(exon_frame[["Name", "Score"]]).drop(like="_b")
    exon_end = PyRanges(sample_end).join(exon_frame[["Name", "Score"]]).drop(like="_b")  

    if exon_start.df.empty and exon_end.df.empty:
        sample_frame['Name'] = sample_frame['Name2'] = sample_frame['Score'] = sample_frame['Score2'] ='None' 
    elif exon_start.df.empty:
        sample_frame = exon_end.df.rename(columns = {'Name':'Name2', 'Score':'Score2'}).filter(items=['ID', 'Name']).drop_duplicates().merge(sample_frame, on=['ID'], how='right')
        sample_frame['Name'] = sample_frame['Score'] = 'None'
    elif exon_end.df.empty: 
        sample_frame = exon_start.df.filter(items=['ID', 'Name', 'Score']).drop_duplicates().merge(sample_frame, on=['ID'], how='right')
        sample_frame['Name2'] = sample_frame['Score2'] = 'None'
    else: 
        sample_frame = exon_start.df.filter(items=['ID', 'Name']).drop_duplicates().merge(sample_frame, on=['ID'], how='right')
        sample_frame = sample_frame.merge(exon_end.df.rename(columns = {'Name':'Name2', 'Score':'Score2'}).filter(items=['ID', 'Name']), on=['ID'], how='left')
    
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
    
    #add column based on overlap with parents
    denovo_start = PyRanges(sample_start).overlap(PyRanges(parent_start))
    denovo_end = PyRanges(sample_end).overlap(PyRanges(parent_end))
    denovo_frame = pd.concat([denovo_start.df, denovo_end.df], axis=0)

    if denovo_frame.empty:
        filtered_sample_frame["Found_in_Parent"] = "False"
        calls = filtered_sample_frame.rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()
    else: 
        denovo_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_frame, on=None, how='left', indicator='Found_in_Parent')
        denovo_filtered_sample_frame['Found_in_Parent'] = np.where(denovo_filtered_sample_frame.Found_in_Parent == 'both', True, False)
        denovo_filtered_sample_frame['CHROM'] = denovo_filtered_sample_frame['CHROM'].map(lambda x: x.lstrip('chr'))
        calls = denovo_filtered_sample_frame.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()

    calls.to_csv(args.outputdirectory+'10xLargeSvInversions'+args.sampleID, sep='\t')

tenxlargesvinversions()
