#!/usr/bin/env python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)
pd.options.display.max_rows = 999



def readsmap(input, args, sv_type):

    colnames = ['SmapEntryID', 'QryContigID', 'RefcontigID1', 'RefcontigID2', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                'RefEndPos', 'Confidence', 'Type', 'XmapID1', 'XmapID2', 'LinkID', 'QryStartIdx', 'QryEndIdx',
                'RefStartIdx', 'RefEndIdx', 'Zygosity', 'Genotype', 'GenotypeGroup', 'RawConfidence',
                'RawConfidenceLeft', 'RawConfidenceRight', 'RawConfidenceCenter', 'SVsize']

    raw_df = pd.read_csv(input, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    raw_df = raw_df[['SmapEntryID', 'RefcontigID1', 'RefcontigID2', 'RefStartPos','RefEndPos', 'QryStartPos', 'QryEndPos', 'Confidence', 'Type', 'Zygosity', 'Genotype']]
    confident_df = raw_df.loc[raw_df['Confidence'] > 0.5] #modulate confidence threshold here
    confident_df  = confident_df[confident_df['Type']==sv_type]
    confident_df['RefcontigID1'] = confident_df['RefcontigID1'].astype(str).str.replace("23", "X")
    confident_df['RefcontigID2'] = confident_df['RefcontigID2'].astype(str).str.replace("24", "Y")

    # calculate SV size
    confident_df['SV_size'] = confident_df['RefEndPos'] - confident_df['RefStartPos'] - confident_df['QryEndPos'] + confident_df['QryStartPos']
    confident_df['SV_size'] = confident_df['SV_size'].abs().round(0)
    confident_df = confident_df.loc[confident_df['SV_size'] >= 1000]

    return(confident_df)





def overlap_length(df):
    overlap = np.maximum(0, np.minimum(df.End, df.End_b) - np.maximum(df.Start, df.Start_b))
    frac = overlap / df.Length
    frac_b = overlap / df.Length_b
    df.insert(df.shape[1], "Overlap", overlap)
    df.insert(df.shape[1], "Fraction", frac)
    df.insert(df.shape[1], "Fraction_b", frac_b)
    return df


def reciprocal_overlap(df1, df2): #takes input PyRanges and refPyRanges
    df1.Length = df1.lengths()
    df2.Length_b = df2.lengths()
    overlap = df1.join(df2)
    overlap = overlap.apply(overlap_length).df
    overlap = overlap.loc[overlap['Fraction'] >= 0.5]
    overlap = overlap.loc[overlap['Fraction_b'] >= 0.5]
    return overlap





def checkRefOverlap(sample_copy, ref_copy, sample_frame):

    overlap_frame = reciprocal_overlap(PyRanges(sample_copy), PyRanges(ref_copy))

    if overlap_frame.empty:
        filtered_sample_frame = sample_frame
    else:
        common = sample_copy.merge(overlap_frame, on=['SmapEntryID'])
        filtered_sample_frame = sample_frame[(~sample_frame.SmapEntryID.isin(common.SmapEntryID))]

    return filtered_sample_frame



def checkParentsOverlap(sample_copy, parent_copy, filtered_sample_frame, args, inheritance):

    if args.type == 'singleton' or (args.type == 'duo' and inheritance == 'Found_in_Father' and args.mother_duo) or  (args.type == 'duo' and inheritance == 'Found_in_Mother' and args.father_duo):
        # Initialize columns and set to -1 if parents file not provided
        filtered_sample_frame[inheritance] = 'None'
        return filtered_sample_frame

    colnames = ['SmapEntryID', 'RefcontigID1', 'RefcontigID2', 'RefStartPos', 'RefEndPos', 'QryStartPos', 'QryEndPos','Confidence', 'Type', 'Zygosity', 'Genotype']

    denovo_parent_frame = reciprocal_overlap(PyRanges(sample_copy), PyRanges(parent_copy))[colnames]

    if denovo_parent_frame.empty:
        filtered_sample_frame[inheritance] = "False"

    else:
        parent_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_parent_frame, on=None, how='left', indicator=inheritance)
        parent_filtered_sample_frame[inheritance] = np.where(parent_filtered_sample_frame[inheritance] == 'both', True,False)
        filtered_sample_frame = parent_filtered_sample_frame.drop_duplicates().reset_index(drop=True)


    return filtered_sample_frame




def exonOverlap(args, df):


    exon_frame = pr.read_bed(args.exons)
    exon_overlap = PyRanges(df).join(exon_frame).drop(like="_b")

    if exon_overlap.df.empty:
        exon_calls = pd.DataFrame()
    else:
        exon_calls = exon_overlap.df.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Name':'gene', 'Score':'OMIM_syndrome'}).drop_duplicates()
        # if args.genelist:
        #     #gene_list = pd.read_csv(args.genelist, sep='\t', names=['Gene'], header=None)
        exon_calls = exon_calls.merge(args.genelist, on=['gene'], how='left')
        exon_calls.fillna(value={'score':0,'normalized_score':0}, inplace=True)
        exon_calls = exon_calls.sort_values(by='score', ascending=False)

    return exon_calls




def BN_deletion(args):

    # Loadsample
    sample_frame = readsmap(args.samplepath, args, 'deletion')

    # Load parent data
    if args.type == 'trio' or args.father_duo:
        father_frame = readsmap(args.fpath, args, 'deletion')
    else:
        father_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1'])

    if args.type == 'trio' or args.mother_duo:
        mother_frame = readsmap(args.mpath, args, 'deletion')
    else:
        mother_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1'])

    # Load reference
    ref_frame = readsmap(args.referencepath, args, 'deletion')


    # Actual work
    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome']  = df.RefStartPos, df.RefEndPos, df['RefcontigID1']

    # Remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)


    # Add column based on overlap with parents
    filtered_sample_frame = checkParentsOverlap(sample_copy, father_copy, filtered_sample_frame, args, 'Found_in_Father')
    filtered_sample_frame = checkParentsOverlap(sample_copy, mother_copy, filtered_sample_frame, args, 'Found_in_Mother')

    # Write the raw output
    filtered_sample_frame.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNano_deletions_raw.txt', sep='\t', index = False)

    # Describe exon overlap
    filtered_sample_frame['Start'], filtered_sample_frame['End'], filtered_sample_frame['Chromosome'] = filtered_sample_frame.RefStartPos, filtered_sample_frame.RefEndPos, filtered_sample_frame['RefcontigID1']
    exon_calls = exonOverlap(args, filtered_sample_frame)

    # Write output files
    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_Bionano_deletions_exons.txt', sep='\t', index = False)

    return filtered_sample_frame, exon_calls




def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-c", "--confidence", help="Give the confidence level cutoff for the sample here", dest="confidence", type=float, default=0.5)
    parser.add_argument("-e", "--exons", help="Give the file with exons intervals, names, and phenotypes here", dest="exons", type=str, required=True)
    parser.add_argument("-g", "--genelist", help="Primary genelist with scores", dest="genelist", type=str)
    parser.add_argument("-t", "--type", help="Specify whether this is a trio, duo, or singleton case", dest="type", type=str)
    parser.add_argument("-F", help="Set this flag if this is a duo case AND only father is sequenced", dest="father_duo", action='store_true')
    parser.add_argument("-M", help="Set this flag if this is a duo case AND only mother is sequenced", dest="mother_duo", action='store_true')
    args = parser.parse_args()

    # Actual function
    BN_deletion(args)

if __name__=="__main__":
    main()
