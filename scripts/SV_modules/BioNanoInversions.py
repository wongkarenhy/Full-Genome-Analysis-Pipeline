#!/usr/bin/env python3.6
  
import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse
pd.set_option('display.max_columns', None)
from .BioNanoTranslocations import checkParentsOverlapTransloInv, geneOverlapTransloInv

def readsmapInv(input):

    colnames = ['SmapEntryID', 'QryContigID', 'RefcontigID1', 'RefcontigID2', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                'RefEndPos', 'Confidence', 'Type', 'XmapID1', 'XmapID2', 'LinkID', 'QryStartIdx', 'QryEndIdx',
                'RefStartIdx', 'RefEndIdx', 'Zygosity', 'Genotype', 'GenotypeGroup', 'RawConfidence',
                'RawConfidenceLeft', 'RawConfidenceRight', 'RawConfidenceCenter', 'SVsize']

    raw_df = pd.read_csv(input, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    raw_df = raw_df[['SmapEntryID', 'RefcontigID1', 'RefcontigID2', 'RefStartPos','RefEndPos', 'QryStartPos', 'QryEndPos', 'Confidence', 'Type', 'Zygosity', 'Genotype']]
    #confident_df = raw_df.loc[raw_df['Confidence'] > 0.5] #modulate confidence threshold here
    confident_df  = raw_df[(raw_df['Type']=='inversion') | (raw_df['Type']=='inversion_paired')]
    confident_df['RefcontigID1'] = confident_df['RefcontigID1'].astype(str).str.replace("23", "X")
    confident_df['RefcontigID2'] = confident_df['RefcontigID2'].astype(str).str.replace("24", "Y")

    return(confident_df)


def BN_inversion(args):

    # Loadsample
    sample_frame = readsmapInv(args.samplepath)

    # Load parent data
    if args.type == 'trio' or args.father_duo:
        father_frame = readsmapInv(args.fpath)
    else:
        father_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1','RefcontigID2'])

    if args.type == 'trio' or args.mother_duo:
        mother_frame = readsmapInv(args.mpath)
    else:
        mother_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1','RefcontigID2'])

    # Load reference
    ref_frame = readsmapInv(args.referencepath)

    # Actual fun
    sample_start, sample_end, mother_start, mother_end, father_start, father_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), mother_frame.copy(), mother_frame.copy(), father_frame.copy(), father_frame.copy(), ref_frame.copy(), ref_frame.copy()

    for df in [sample_start, mother_start, father_start,
               ref_start]:  # create an interval for the translocation start point
        df['Start'], df['End'], df['Chromosome'] = df.RefStartPos - 20000, df.RefStartPos + 20000, df['RefcontigID1']

    for df in [sample_end, mother_end, father_end, ref_end]:  # create an interval for the translocation end point
        df['Start'], df['End'], df['Chromosome'] = df.RefEndPos - 20000, df.RefEndPos + 20000, df['RefcontigID2']


    # Overlap start and end points with genes separately
    sample_frame = geneOverlapTransloInv(args, sample_start, sample_end, sample_frame)

    # Remove anything that overlaps with the reference
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

    # Add column based on overlap with parent
    filtered_sample_frame = checkParentsOverlapTransloInv(filtered_sample_frame, sample_start, father_start, sample_end, father_end, args, 'Found_in_Father')
    filtered_sample_frame = checkParentsOverlapTransloInv(filtered_sample_frame, sample_start, mother_start, sample_end, mother_end, args, 'Found_in_Mother')

    # Clean up dataframe
    if 'Chromosome' in list(filtered_sample_frame.columns):
        calls = filtered_sample_frame.drop(columns=['Chromosome', 'Start', 'End']).rename(columns={'Name': 'Gene', 'Name2': 'Gene2', 'Score': 'OMIM_syndrome', 'Score2': 'OMIM_syndrome2'}).drop_duplicates()
    else:
        calls = filtered_sample_frame.rename(columns={'Name': 'Gene', 'Name2': 'Gene2', 'Score': 'OMIM_syndrome', 'Score2': 'OMIM_syndrome2'}).drop_duplicates()

    # Write output
    calls.to_csv(args.outputdirectory + '/confident_set/' + args.sampleID + '_Bionano_inversions.txt', sep='\t', index = False)



def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-c", "--confidence", help="Give the confidence level cutoff for the sample here", dest="confidence", type=str, default=0.5)
    parser.add_argument("-e", "--genes", help="Give the file with genes intervals, names, and phenotypes here", dest="genes", type=str, required=True)
    parser.add_argument("-t", "--type", help="Specify whether this is a trio, duo, or singleton case", dest="type", type=str)
    parser.add_argument("-F", help="Set this flag if this is a duo case AND only father is sequenced", dest="father_duo", action='store_true')
    parser.add_argument("-M", help="Set this flag if this is a duo case AND only mother is sequenced", dest="mother_duo", action='store_true')
    args = parser.parse_args()


    # Actual function
    BN_inversion(args)


if __name__=="__main__":
    main()
