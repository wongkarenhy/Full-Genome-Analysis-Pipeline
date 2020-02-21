#!/usr/bin/env python3.6
  
import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse
pd.set_option('display.max_columns', None)
from .BioNanoTranslocations import checkParentsOverlapTransloInv

def readsmapInv(input):

    colnames = ['SmapEntryID', 'QryContigID', 'RefcontigID1', 'RefcontigID2', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                'RefEndPos', 'Confidence', 'Type', 'XmapID1', 'XmapID2', 'LinkID', 'QryStartIdx', 'QryEndIdx',
                'RefStartIdx', 'RefEndIdx', 'Zygosity', 'Genotype', 'GenotypeGroup', 'RawConfidence',
                'RawConfidenceLeft', 'RawConfidenceRight', 'RawConfidenceCenter', 'SVsize']

    raw_df = pd.read_csv(input, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    raw_df = raw_df[['SmapEntryID', 'RefcontigID1', 'RefcontigID2', 'RefStartPos','RefEndPos', 'QryStartPos', 'QryEndPos', 'Confidence', 'Type', 'Zygosity', 'Genotype']]
    #confident_df = raw_df.loc[raw_df['Confidence'] > 0.5] #modulate confidence threshold here
    confident_df  = raw_df[raw_df['Type']=='inversion']

    return(confident_df)


def BN_inversion(args):

    #loadsample
    sample_frame = readsmapInv(args.samplepath)

    #loadparent
    father_frame = readsmapInv(args.fpath)
    mother_frame = readsmapInv(args.mpath)

    #load reference
    ref_frame = readsmapInv(args.referencepath)

    # Actual fun
    sample_start, sample_end, mother_start, mother_end, father_start, father_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), mother_frame.copy(), mother_frame.copy(), father_frame.copy(), father_frame.copy(), ref_frame.copy(), ref_frame.copy()

    for df in [sample_start, mother_start, father_start,
               ref_start]:  # create an interval for the translocation start point
        df['Start'], df['End'], df['Chromosome'] = df.RefStartPos - 100000, df.RefStartPos + 100000, df['RefcontigID1']

    for df in [sample_end, mother_end, father_end, ref_end]:  # create an interval for the translocation end point
        df['Start'], df['End'], df['Chromosome'] = df.RefEndPos - 100000, df.RefEndPos + 100000, df['RefcontigID2']


    #overlap start and end points with exons separately  
    exon_frame = pr.read_bed(args.exons)
    exon_start = PyRanges(sample_start).join(exon_frame[["Name", "Score"]]).drop(like="_b")
    exon_end = PyRanges(sample_end).join(exon_frame[["Name", "Score"]]).drop(like="_b")

    if exon_start.df.empty and exon_end.df.empty:
        sample_frame['Name'] = sample_frame['Name2'] = sample_frame['Score'] = sample_frame['Score2'] ='None'
    elif exon_start.df.empty:
        sample_frame = exon_end.df.rename(columns = {'Name':'Name2', 'Score':'Score2'}).filter(items=['SmapEntryID', 'Name']).drop_duplicates().merge(sample_frame, on=['SmapEntryID'], how='right')
        sample_frame['Name'] = sample_frame['Score'] = 'None'
    elif exon_end.df.empty:
        sample_frame = exon_start.df.filter(items=['SmapEntryID', 'Name', 'Score']).drop_duplicates().merge(sample_frame, on=['SmapEntryID'], how='right')
        sample_frame['Name2'] = sample_frame['Score2'] = 'None'
    else:
        sample_frame = exon_start.df.filter(items=['SmapEntryID', 'Name', 'Score']).drop_duplicates().merge(sample_frame, on=['SmapEntryID'], how='right')
        sample_frame = sample_frame.merge(exon_end.df.rename(columns = {'Name':'Name2', 'Score':'Score2'}).filter(items=['SmapEntryID', 'Name2', 'Score2']), on=['SmapEntryID'], how='left')

    #remove anything that overlaps with the reference
    overlap_start = PyRanges(sample_start).overlap(PyRanges(ref_start))
    overlap_end = PyRanges(sample_end).overlap(PyRanges(ref_end))
    if overlap_start.df.empty and overlap_end.df.empty:
        filtered_sample_frame = sample_frame
    else:
        overlap_frame = overlap_start.df.merge(overlap_end.df, on=['SmapEntryID'])
        if overlap_frame.empty:
            filtered_sample_frame = sample_frame
        else:
            common = sample_frame.merge(overlap_frame.filter(items=['SmapEntryID']),on=['SmapEntryID'])
            filtered_sample_frame = sample_frame[(~sample_frame.SmapEntryID.isin(common.SmapEntryID))]

    #add column based on overlap with parents
    calls = checkParentsOverlapTransloInv(filtered_sample_frame, sample_start, father_start, mother_start, sample_end, father_end, mother_end)

    # Write output
    calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_Bionano_inversions.txt', sep='\t', index = False)



def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-c", "--confidence", help="Give the confidence level cutoff for the sample here", dest="confidence", type=str, default=0.5)
    parser.add_argument("-e", "--exons", help="Give the file with exons intervals, names, and phenotypes here", dest="exons", type=str, required=True)
    args = parser.parse_args()


    # Actual function
    BN_inversion(args)


if __name__=="__main__":
    main()
