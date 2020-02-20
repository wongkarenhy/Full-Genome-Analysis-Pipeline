#!/usr/bin/env python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse


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



def readsmap(input, args, sv_type):

    colnames = ['SmapEntryID', 'QryContigID', 'RefcontigID1', 'RefcontigID2', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                'RefEndPos', 'Confidence', 'Type', 'XmapID1', 'XmapID2', 'LinkID', 'QryStartIdx', 'QryEndIdx',
                'RefStartIdx', 'RefEndIdx', 'Zygosity', 'Genotype', 'GenotypeGroup', 'RawConfidence',
                'RawConfidenceLeft', 'RawConfidenceRight', 'RawConfidenceCenter', 'SVsize']

    raw_df = pd.read_csv(input, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    raw_df = raw_df[['SmapEntryID', 'RefcontigID1', 'RefcontigID2', 'RefStartPos','RefEndPos', 'QryStartPos', 'QryEndPos', 'Confidence', 'Type', 'Zygosity', 'Genotype']]
    confident_df = raw_df.loc[raw_df['Confidence'] > 0.5] #modulate confidence threshold here
    confident_df  = confident_df[confident_df['Type']==sv_type]

    return(confident_df)



def cytobandOverlap(args, df):

    #describe cytoband overlap
    cytoband_frame = pr.read_bed(args.cytobands)
    cytoband_overlap = PyRanges(df).join(cytoband_frame).drop(like="_b")
    cytoband_calls = cytoband_overlap.df
    if not cytoband_calls.empty:
        cytoband_calls["Cytoband"] = cytoband_calls[['Chromosome', 'Name']].apply(lambda x: ''.join(x), axis=1)
        cytoband_calls = cytoband_calls.drop(columns = ['Chromosome', 'Start', 'End', 'Name'])

    return cytoband_calls



def checkRefOverlap(sample_copy, ref_copy, sample_frame):

    overlap_frame = reciprocal_overlap(PyRanges(sample_copy), PyRanges(ref_copy))

    if overlap_frame.empty:
        filtered_sample_frame = sample_frame
    else:
        common = sample_copy.merge(overlap_frame, on=['SmapEntryID'])
        filtered_sample_frame = sample_frame[(~sample_frame.SmapEntryID.isin(common.SmapEntryID))]

    return filtered_sample_frame



def checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame):

    denovo_f_frame, denovo_m_frame = reciprocal_overlap(PyRanges(sample_copy), PyRanges(father_copy)), reciprocal_overlap(PyRanges(sample_copy), PyRanges(mother_copy))

    if denovo_f_frame.empty:
        filtered_sample_frame["Found_in_Father"] = "False"
    else:
        filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_f_frame, on=None, how='left', indicator='Found_in_Father')
        filtered_sample_frame['Found_in_Father'] = np.where(filtered_sample_frame.Found_in_Father == 'both', True, False)
        filtered_sample_frame = filtered_sample_frame.drop(columns=['Start', 'End', 'Chromosome', 'Length', 'SmapEntryID_b', 'RefcontigID1_b', 'RefcontigID2_b', 'RefStartPos_b', 'RefEndPos_b', 'QryStartPos_b', 'QryEndPos_b', 'Confidence_b', 'Type_b', 'Zygosity_b', 'Genotype_b', 'Start_b', 'End_b', 'Length_b', 'Overlap', 'Fraction', 'Fraction_b']).drop_duplicates()


    if denovo_m_frame.empty:
        filtered_sample_frame["Found_in_Mother"] = "False"
        denovo_filtered_sample_frame = filtered_sample_frame

    else:
        denovo_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_m_frame, on=None, how='left', indicator='Found_in_Mother')
        denovo_filtered_sample_frame['Found_in_Mother'] = np.where(denovo_filtered_sample_frame.Found_in_Mother == 'both', True, False)
        denovo_filtered_sample_frame = denovo_filtered_sample_frame.drop(columns=['Start', 'End', 'Chromosome', 'Length', 'SmapEntryID_b', 'RefcontigID1_b', 'RefcontigID2_b', 'RefStartPos_b', 'RefEndPos_b', 'QryStartPos_b', 'QryEndPos_b', 'Confidence_b', 'Type_b', 'Zygosity_b', 'Genotype_b', 'Start_b', 'End_b', 'Length_b', 'Overlap', 'Fraction', 'Fraction_b']).drop_duplicates()

    df = denovo_filtered_sample_frame

    # Make proper columns for the next step (exon overlap)
    df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']


    return df




def exonOverlap(args, df):


    exon_frame = pr.read_bed(args.exons)
    exon_overlap = PyRanges(df).join(exon_frame).drop(like="_b")
    if exon_overlap.df.empty:
        exon_calls = pd.DataFrame()
    else:
        exon_calls = exon_overlap.df.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Name':'gene', 'Score':'Phenotype'}).drop_duplicates()

        # if args.genelist:
        #     #gene_list = pd.read_csv(args.genelist, sep='\t', names=['Gene'], header=None)
        exon_calls = exon_calls.merge(args.genelist, on=['gene'])
        exon_calls = exon_calls.sort_values(by='score', ascending=False)

    return exon_calls




def BN_deletion(args):

    #loadsample
    sample_frame = readsmap(args.samplepath, args, 'deletion')

    #loadparent
    mother_frame = readsmap(args.mpath, args, 'deletion')
    father_frame = readsmap(args.fpath, args, 'deletion')

    # #load reference
    ref_frame = readsmap(args.referencepath, args, 'deletion')


    # Actual work
    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome']  = df.RefStartPos, df.RefEndPos, df['RefcontigID1']

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)


    #add column based on overlap with parents
    df = checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame)

    #cytoband_calls = cytobandOverlap(args, df)
    #df.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNanoDeletions_cytobands.txt', sep='\t', index = False)

    #describe exon overlap
    exon_calls = exonOverlap(args, df)

    # Write output files
    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNanoDeletions_exons.txt', sep='\t', index = False)

    return df




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
    args = parser.parse_args()

    # Actual function
    BN_deletion(args)

if __name__=="__main__":
    main()
