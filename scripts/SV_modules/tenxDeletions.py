#!/usr/bin/env python3.6

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import allel
import argparse
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)
pd.options.display.max_rows = 999



def overlap_length(df):

    overlap = np.maximum(0, np.minimum(df.End, df.End_b) - np.maximum(df.Start, df.Start_b))
    frac = overlap / df.Length
    frac_b = overlap / df.Length_b
    df.insert(df.shape[1], "Overlap", overlap)
    df.insert(df.shape[1], "Fraction", frac)
    df.insert(df.shape[1], "Fraction_b", frac_b)

    return df


def reciprocal_overlap(overlap): #takes input PyRanges and refPyRanges

    # df1.Length = df1.lengths()
    # df2.Length_b = df2.lengths()
    # overlap = df1.join(df2)
    overlap = overlap.apply(overlap_length).df
    overlap = overlap.loc[overlap['Fraction'] >= 0.5]
    overlap = overlap.loc[overlap['Fraction_b'] >= 0.5]
    return overlap



def read10x(input):

    df = allel.vcf_to_dataframe(input, fields=['variants/CHROM', 'variants/POS', 'variants/ID', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER_PASS', 'variants/END', 'variants/SVLEN'])
    # scores_cutoff = np.mean(df.QUAL) - 2*np.std(df.QUAL)
    # df = df.loc[df['FILTER_PASS']==True]
    # df = df.loc[df['QUAL']>scores_cutoff]
    df.reset_index(inplace = True, drop=True)
    df['CHROM'] = df['CHROM'].map(lambda x: x.lstrip('chr'))

    return(df)




def checkRefOverlap(sample_copy, ref_copy, sample_frame):


    sample_copy_py, ref_copy_py = PyRanges(sample_copy), PyRanges(ref_copy)
    sample_copy_py.Length, ref_copy_py.Length_b = sample_copy_py.lengths(), ref_copy_py.lengths()
    overlap = sample_copy_py.join(ref_copy_py)

    #overlap = PyRanges(sample_copy).join(PyRanges(ref_copy))
    #print(overlap.head())
    overlap_frame = reciprocal_overlap(overlap)

    if overlap_frame.empty:
        filtered_sample_frame = sample_frame
    else:
        common = sample_copy.merge(overlap_frame, on=['ID'])
        filtered_sample_frame = sample_frame[(~sample_frame.ID.isin(common.ID))]

    return filtered_sample_frame




def checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame):

    colnames = ['CHROM', 'POS', 'ID', 'REF', "ALT_1", 'ALT_2', "ALT_3", "QUAL", 'FILTER_PASS', "END", "SVLEN"]

    sample_copy_py, father_copy_py, mother_copy_py = PyRanges(sample_copy), PyRanges(father_copy), PyRanges(mother_copy)
    sample_copy_py.Length, father_copy_py.Length_b, mother_copy_py.Length_b = sample_copy_py.lengths(), father_copy_py.lengths(), mother_copy_py.lengths()

    overlap_sample_father = sample_copy_py.join(father_copy_py)
    overlap_sample_mother = sample_copy_py.join(mother_copy_py)

    #denovo_f_frame, denovo_m_frame = reciprocal_overlap(PyRanges(sample_copy), PyRanges(father_copy))[colnames], reciprocal_overlap(PyRanges(sample_copy), PyRanges(mother_copy))[colnames]

    denovo_f_frame = reciprocal_overlap(overlap_sample_father)[colnames]
    denovo_m_frame = reciprocal_overlap(overlap_sample_mother)[colnames]

    if denovo_f_frame.empty:
        filtered_sample_frame["Found_in_Father"] = "False"
    if denovo_m_frame.empty:
        filtered_sample_frame["Found_in_Mother"] = "False"

    elif not (denovo_f_frame.empty and denovo_m_frame.empty):
        f_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_f_frame, on=None, how='left',
                                           indicator='Found_in_Father')
        f_filtered_sample_frame['Found_in_Father'] = np.where(f_filtered_sample_frame.Found_in_Father == 'both', True, False)
        f_filtered_sample_frame = f_filtered_sample_frame.drop_duplicates().reset_index(drop=True)

        m_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_m_frame, on=None, how='left',
                                           indicator='Found_in_Mother')
        m_filtered_sample_frame['Found_in_Mother'] = np.where(m_filtered_sample_frame.Found_in_Mother == 'both', True,False)
        m_filtered_sample_frame = m_filtered_sample_frame.drop_duplicates().reset_index(drop=True)
        f_filtered_sample_frame['Found_in_Mother'] = m_filtered_sample_frame['Found_in_Mother']
        filtered_sample_frame = f_filtered_sample_frame

    # df = filtered_sample_frame
    # df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']


    # denovo_f_frame = reciprocal_overlap(overlap_sample_father)
    # denovo_m_frame = reciprocal_overlap(overlap_sample_mother)
    #
    # if denovo_f_frame.empty:
    #     filtered_sample_frame["Found_in_Father"] = "False"
    # else:
    #     filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_f_frame, on=None, how='left', indicator='Found_in_Father')
    #     filtered_sample_frame['Found_in_Father'] = np.where(filtered_sample_frame.Found_in_Father == 'both', True, False)
    #     filtered_sample_frame = filtered_sample_frame.drop(columns=['Start', 'End', 'Chromosome', 'Length', 'Length_b', 'Overlap', 'Fraction', 'Fraction_b', 'CHROM_b', 'POS_b', 'ID_b', 'REF_b', 'ALT_1_b','ALT_2_b','ALT_3_b', 'QUAL_b','FILTER_PASS_b', 'END_b','SVLEN_b','Start_b', 'End_b']).drop_duplicates()
    #
    #
    # if denovo_m_frame.empty:
    #     filtered_sample_frame["Found_in_Mother"] = "False"
    #     denovo_filtered_sample_frame = filtered_sample_frame
    #
    # else:
    #     denovo_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_m_frame, on=None, how='left', indicator='Found_in_Mother')
    #     denovo_filtered_sample_frame['Found_in_Mother'] = np.where(denovo_filtered_sample_frame.Found_in_Mother == 'both', True, False)
    #     denovo_filtered_sample_frame = denovo_filtered_sample_frame.drop(columns=['Start', 'End', 'Chromosome', 'Length', 'Length_b', 'Overlap', 'Fraction', 'Fraction_b', 'CHROM_b', 'POS_b', 'ID_b', 'REF_b', 'ALT_1_b','ALT_2_b','ALT_3_b', 'QUAL_b','FILTER_PASS_b', 'END_b','SVLEN_b','Start_b', 'End_b']).drop_duplicates()
    #
    # df = denovo_filtered_sample_frame
    #
    # # Make proper columns for the next step (exon overlap)
    # df['CHROM'] = df['CHROM'].map(lambda x: x.lstrip('chr'))
    # df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']

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




def tenxdeletions(args):

    #load sample data
    sample_frame = read10x(args.samplepath)

    if not args.singleton:
        # load parent data
        father_frame = read10x(args.fpath)
        mother_frame = read10x(args.mpath)
    else:
        mother_frame = pd.DataFrame(columns=['POS', 'END', 'CHROM'])
        father_frame = pd.DataFrame(columns=['POS', 'END', 'CHROM'])


    #load reference data
    ref_frame = read10x(args.referencepath)

    sample_copy, mother_copy, father_copy, ref_copy = sample_frame.copy(), mother_frame.copy(), father_frame.copy(), ref_frame.copy()

    for df in [sample_copy, mother_copy, father_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']

    #remove anything that overlaps with the reference
    filtered_sample_frame = checkRefOverlap(sample_copy, ref_copy, sample_frame)

    #add column based on overlap with parents
    if not args.singleton:
        df = checkParentsOverlap(sample_copy, father_copy, mother_copy, filtered_sample_frame)
    else:
        df = filtered_sample_frame

    #describe exon overlap
    df['Start'], df['End'], df['Chromosome'] = df.POS, df.END, df['CHROM']
    exon_calls = exonOverlap(args, df)

    # Write final output
    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_10x_deletions_exons.txt', sep='\t', index = False)

    return df, exon_calls


def main():

    parser = argparse.ArgumentParser(description='Sample and Path arguments.')
    parser.add_argument("-i", "--sampleID", help="Give the sample ID", dest="sampleID", type=str, required=True)
    parser.add_argument("-s", "--samplepath", help="Give the full path to the sample file", dest="samplepath", type=str, required=True)
    parser.add_argument("-f", "--fpath", help="Give the full path to the father's file", dest="fpath", type=str, required=True)
    parser.add_argument("-m", "--mpath", help="Give the full path to the mother's file", dest="mpath", type=str, required=True)
    parser.add_argument("-r", "--referencepath", help="Give the full path to the reference file", dest="referencepath", type=str, required=True)
    parser.add_argument("-o", "--outputdirectory", help="Give the directory path for the output file", dest="outputdirectory", type=str, required=True)
    parser.add_argument("-e", "--exons", help="Give the file with exons intervals, names, and phenotypes here", dest="exons", type=str, required=True)
    parser.add_argument("-g", "--genelist", help="Primary genelist with scores", dest="genelist", type=str)
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()

    # Actual function
    tenxdeletions(args)

if __name__=="__main__":
    main()

