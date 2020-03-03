#!/usr/bin/env python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse
from .BioNanoDeletions import reciprocal_overlap, checkRefOverlap, checkParentsOverlap, exonOverlap, overlap_length


def readsmapDup(input):

    colnames = ['SmapEntryID', 'QryContigID', 'RefcontigID1', 'RefcontigID2', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                'RefEndPos', 'Confidence', 'Type', 'XmapID1', 'XmapID2', 'LinkID', 'QryStartIdx', 'QryEndIdx',
                'RefStartIdx', 'RefEndIdx', 'Zygosity', 'Genotype', 'GenotypeGroup', 'RawConfidence',
                'RawConfidenceLeft', 'RawConfidenceRight', 'RawConfidenceCenter', 'SVsize']

    raw_df = pd.read_csv(input, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    raw_df = raw_df[['SmapEntryID', 'RefcontigID1', 'RefcontigID2', 'RefStartPos','RefEndPos', 'QryStartPos', 'QryEndPos', 'Confidence', 'Type', 'Zygosity', 'Genotype']]
    #confident_df = raw_df.loc[raw_df['Confidence'] > 0.5] #modulate confidence threshold here
    confident_df  = raw_df[raw_df['Type'].str.contains('duplication')]
    confident_df['RefcontigID1'] = confident_df['RefcontigID1'].astype(str).str.replace("23", "X")
    confident_df['RefcontigID2'] = confident_df['RefcontigID2'].astype(str).str.replace("24", "Y")


    # calculate SV size
    confident_df['SV_size'] = confident_df['RefEndPos'] - confident_df['RefStartPos'] - confident_df['QryEndPos'] + confident_df['QryStartPos']
    confident_df['SV_size'] = confident_df['SV_size'].abs().round(0)
    confident_df = confident_df.loc[confident_df['SV_size'] >= 1000]

    return(confident_df)



def BN_duplication(args):

    #loadsample
    sample_frame = readsmapDup(args.samplepath)

    # Some old BN pipeline doesn't call duplication
    if sample_frame.empty:
        return None, None

    #loadparent
    father_frame = readsmapDup(args.fpath)
    mother_frame = readsmapDup(args.mpath)

    #load reference
    ref_frame = readsmapDup(args.referencepath)

    # Actual fun
    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome']  = df.RefStartPos, df.RefEndPos, df['RefcontigID1']

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)

    #add column based on overlap with parents
    if not args.singleton:
        df = checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame)
    else:
        df = filtered_sample_frame

    #describe exon overlap
    df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']
    exon_calls = exonOverlap(args, df)


    #df.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNanoDuplications_cytobands.txt', sep='\t', index = False)
    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_Bionano_duplications_exons.txt', sep='\t', index = False)

    return df, exon_calls


                
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
    parser.add_argument("-g", "--genelist", help="Primary genelist with scores", dest="genelist", type=str)
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()


    # Actual function
    BN_duplication(args)


if __name__=="__main__":
    main()

