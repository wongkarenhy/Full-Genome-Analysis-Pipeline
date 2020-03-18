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



def checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame):
    colnames = ['SmapEntryID', 'RefcontigID1', 'RefcontigID2', 'RefStartPos', 'RefEndPos', 'QryStartPos', 'QryEndPos','Confidence', 'Type', 'Zygosity', 'Genotype']

    denovo_f_frame, denovo_m_frame = reciprocal_overlap(PyRanges(sample_copy), PyRanges(father_copy))[colnames], \
                                     reciprocal_overlap(PyRanges(sample_copy), PyRanges(mother_copy))[colnames]

    #print(denovo_f_frame)
    if denovo_f_frame.empty:
        filtered_sample_frame["Found_in_Father"] = "False"
    if denovo_m_frame.empty:
        filtered_sample_frame["Found_in_Mother"] = "False"

    elif not (denovo_f_frame.empty and denovo_m_frame.empty):
        f_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_f_frame, on=None, how='left',
                                           indicator='Found_in_Father')
        f_filtered_sample_frame['Found_in_Father'] = np.where(f_filtered_sample_frame.Found_in_Father == 'both', True,
                                                              False)
        f_filtered_sample_frame = f_filtered_sample_frame.drop_duplicates().reset_index(drop=True)

        m_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_m_frame, on=None, how='left',
                                           indicator='Found_in_Mother')
        m_filtered_sample_frame['Found_in_Mother'] = np.where(m_filtered_sample_frame.Found_in_Mother == 'both', True,
                                                              False)
        m_filtered_sample_frame = m_filtered_sample_frame.drop_duplicates().reset_index(drop=True)
        f_filtered_sample_frame['Found_in_Mother'] = m_filtered_sample_frame['Found_in_Mother']
        filtered_sample_frame = f_filtered_sample_frame

    # df = filtered_sample_frame
    # df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']

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
        exon_calls = exon_calls.merge(args.genelist, on=['gene'])
        exon_calls = exon_calls.sort_values(by='score', ascending=False)

    return exon_calls




def BN_deletion(args):

    #loadsample
    sample_frame = readsmap(args.samplepath, args, 'deletion')

    if not args.singleton:
        #loadparent
        mother_frame = readsmap(args.mpath, args, 'deletion')
        father_frame = readsmap(args.fpath, args, 'deletion')
    else:
        mother_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1'])
        father_frame = pd.DataFrame(columns=['RefStartPos', 'RefEndPos', 'RefcontigID1'])


    # #load reference
    ref_frame = readsmap(args.referencepath, args, 'deletion')


    # Actual work
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
    df.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNano_deletions_raw.txt', sep='\t', index = False)
    df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']
    exon_calls = exonOverlap(args, df)

    # Write output files
    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_Bionano_deletions_exons.txt', sep='\t', index = False)

    return df, exon_calls




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
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()

    # Actual function
    BN_deletion(args)

if __name__=="__main__":
    main()
