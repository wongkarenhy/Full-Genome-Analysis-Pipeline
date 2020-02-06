#!/usr/local/bin/python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import allel
import argparse

def tenxlargesvduplications():

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
    sample_frame = sample_frame.loc[sample_frame['FILTER_PASS']==True].loc[sample_frame['QUAL']>scores_cutoff].loc[sample_frame['ALT_1'].str.contains("DUP")]

    #load parent data 
    father_frame = allel.vcf_to_dataframe(args.fpath, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END'])
    scores_cutoff = np.mean(father_frame.QUAL) - 2*np.std(father_frame.QUAL)
    father_frame = father_frame.loc[father_frame['FILTER_PASS']==True].loc[father_frame['QUAL']>scores_cutoff].loc[father_frame['ALT_1'].str.contains("DUP")]
    mother_frame = allel.vcf_to_dataframe(args.mpath, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END'])
    scores_cutoff = np.mean(mother_frame.QUAL) - 2*np.std(mother_frame.QUAL)
    mother_frame = mother_frame.loc[mother_frame['FILTER_PASS']==True].loc[mother_frame['QUAL']>scores_cutoff].loc[mother_frame['ALT_1'].str.contains("DUP")]
    parent_frame = pd.concat([mother_frame, father_frame])

    #load reference data
    ref_frame = allel.vcf_to_dataframe(args.referencepath, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END'])
    scores_cutoff = np.mean(ref_frame.QUAL) - 2*np.std(ref_frame.QUAL)
    ref_frame = ref_frame.loc[ref_frame['FILTER_PASS']==True]
    ref_frame = ref_frame.loc[ref_frame['QUAL']>scores_cutoff]
    ref_frame = ref_frame.loc[ref_frame['ALT_1'].str.contains("DUP")]

    sample_copy, parent_copy, ref_copy = sample_frame.copy(), parent_frame.copy(), ref_frame.copy()

    for df in [sample_copy, parent_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']

    #remove anything that overlaps with the reference
    overlap_frame = PyRanges(sample_copy).overlap(PyRanges(ref_copy))

    if overlap_frame.df.empty:
        filtered_sample_frame = sample_frame
    else:
        common = sample_copy.merge(overlap_frame.df, on=['ID'])
        filtered_sample_frame = sample_frame[(~sample_frame.ID.isin(common.ID))]

    #add column based on overlap with parents 
    denovo_frame = PyRanges(sample_copy).overlap(PyRanges(parent_copy))

    if denovo_frame.df.empty:
        filtered_sample_frame["Found_in_Parent"] = "False"
        denovo_filtered_sample_frame = filtered_sample_frame

    else:
        denovo_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_frame.df, on=None, how='left', indicator='Found_in_Parent')
        denovo_filtered_sample_frame['Found_in_Parent'] = np.where(denovo_filtered_sample_frame.Found_in_Parent == 'both', True, False)
        denovo_filtered_sample_frame = denovo_filtered_sample_frame.drop_duplicates()

    df = denovo_filtered_sample_frame
    df['CHROM'] = df['CHROM'].map(lambda x: x.lstrip('chr'))
    df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']
    

    #describe cytoband overlap
    cytoband_frame = pr.read_bed(args.cytobands)
    cytoband_overlap = PyRanges(df).join(cytoband_frame).drop(like="_b")
    cytoband_calls = cytoband_overlap.df
    if not cytoband_calls.empty: 
        cytoband_calls["Cytoband"] = cytoband_calls[['Chromosome', 'Name']].apply(lambda x: ''.join(x), axis=1)
        cytoband_calls = cytoband_calls.drop(columns = ['Chromosome', 'Start', 'End', 'Name'])


    #describe exon overlap
    exon_frame = pr.read_bed(args.exons)
    exon_overlap = PyRanges(df).join(exon_frame).drop(like="_b")
    if exon_overlap.df.empty:
        exon_calls = pd.DataFrame()
    else:
        exon_calls = exon_overlap.df.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Name':'Gene', 'Score':'Phenotype'}).drop_duplicates()

    cytoband_calls.to_csv(args.outputdirectory+'10xLargeSvDuplications'+args.sampleID+'cytobands', sep='\t')
    exon_calls.to_csv(args.outputdirectory+'10xLargeSvDuplications'+args.sampleID+'exons', sep='\t')

tenxlargesvduplications()
