#!/usr/bin/env python3

import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse
from BioNanoDeletions import readsmap, cytobandOverlap


def duplication(args):


    # colnames=['SmapEntryID','QryContigID','RefcontigID1','RefcontigID2','QryStartPos','QryEndPos','RefStartPos','RefEndPos','Confidence','Type','XmapID1','XmapID2','LinkID','QryStartIdx','QryEndIdx','RefStartIdx','RefEndIdx','Zygosity','Genotype','GenotypeGroup','RawConfidence','RawConfidenceLeft','RawConfidenceRight','RawConfidenceCenter', 'SVsize']
    #
    # #loadsample
    sample_frame = readsmap(args.samplepath, args, 'duplication')
    # raw_sf = pd.read_csv(args.samplepath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # confident_sf = raw_sf.loc[raw_sf['Confidence'] > args.confidence] #modulate confidence threshold here
    # sample_frame  = confident_sf[confident_sf['Type']=='duplication']
    #
    # #loadparent
    father_frame = readsmap(args.fpath, args, 'duplication')
    mother_frame = readsmap(args.mpath, args, 'duplication')
    parent_frame = pd.concat([mother_frame, father_frame])
    # mother_pf = pd.read_csv(args.mpath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # father_pf = pd.read_csv(args.fpath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # raw_pf = pd.concat([mother_pf, father_pf])
    # confident_pf = raw_pf.loc[raw_pf['Confidence'] > 0.3]
    # parent_frame = confident_pf[confident_pf['Type']=='duplication']
    #
    # #load reference
    ref_frame = readsmap(args.referencepath, args, 'duplication')
    # raw_rf = pd.read_csv(args.referencepath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # confident_rf = raw_rf.loc[raw_rf['Confidence'] > 0.3]
    # del_confident_rf = confident_rf[confident_rf['Type']=='duplication']
    # ref_frame = del_confident_rf
    
    sample_copy, parent_copy, ref_copy = sample_frame.copy(), parent_frame.copy(), ref_frame.copy()

    for df in [sample_copy, parent_copy, ref_copy]:
        df['Start'], df['End'], df['Chromosome']  = df.RefStartPos, df.RefEndPos, df['RefcontigID1']

    #remove anything that overlaps with the reference
    overlap_frame = PyRanges(sample_copy).overlap(PyRanges(ref_copy))

    if overlap_frame.df.empty:
        filtered_sample_frame = sample_frame
    else:
        common = sample_copy.merge(overlap_frame.df, on=['SmapEntryID'])
        filtered_sample_frame = sample_frame[(~sample_frame.SmapEntryID.isin(common.SmapEntryID))]

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
    df['Start'], df['End'], df['Chromosome'] = df.RefStartPos, df.RefEndPos, df['RefcontigID1']

    # #describe cytoband overlap
    cytoband_frame = pr.read_bed(args.cytobands)
    cytoband_overlap = PyRanges(df).join(cytoband_frame).drop(like="_b")
    cytoband_calls = cytoband_overlap.df
    cytoband_calls["Cytoband"] = cytoband_calls[['Chromosome', 'Name']].apply(lambda x: ''.join(x), axis=1)
    cytoband_calls = cytoband_calls.drop(columns = ['Chromosome', 'Start', 'End', 'Name'])
    #cytoband_calls = cytobandOverlap(args, df)

    #describe exon overlap
    exon_frame = pr.read_bed(args.exons)
    exon_overlap = PyRanges(df).join(exon_frame).drop(like="_b")
    if exon_overlap.df.empty:
        exon_calls = pd.Dataframe()
    else:
        exon_calls = exon_overlap.df.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Name':'Gene', 'Score':'Phenotype'}).drop_duplicates()


    cytoband_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNanoDuplications_cytobands.txt', sep='\t', index = False)
    exon_calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNanoDuplications_exons.txt', sep='\t', index = False)
    #print(exon_calls)

                
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
    parser.add_argument("-y", "--cytobands", help="Give the BED file with cytoband intervals", dest="cytobands", type=str, required=True)
    parser.add_argument("-g", "--genelist", help="Primary genelist", dest="genelist", type=str)
    args = parser.parse_args()


    # Actual function
    duplication(args)


if __name__=="__main__":
    main()

