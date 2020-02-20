#!/usr/bin/env python3
  
import pandas as pd
import os
import pyranges as pr
from pyranges import PyRanges
from io import StringIO
import numpy as np
import argparse
pd.set_option('display.max_columns', None)
from BioNanoDeletions import readsmap

def inversion(args):


    # colnames=['SmapEntryID','QryContigID','RefcontigID1','RefcontigID2','QryStartPos','QryEndPos','RefStartPos','RefEndPos','Confidence','Type','XmapID1','XmapID2','LinkID','QryStartIdx','QryEndIdx','RefStartIdx','RefEndIdx','Zygosity','Genotype','GenotypeGroup','RawConfidence','RawConfidenceLeft','RawConfidenceRight','RawConfidenceCenter', 'SVsize']
    #
    #loadsample
    sample_frame = readsmap(args.samplepath, args, 'inversion')
    # raw_sf = pd.read_csv(args.samplepath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # confident_sf = raw_sf.loc[raw_sf['Confidence'] > float(args.confidence)] #modulate confidence threshold here
    # sample_frame  = confident_sf[confident_sf['Type']=='inversion']
    #
    #loadparent
    father_frame = readsmap(args.fpath, args, 'inversion')
    mother_frame = readsmap(args.mpath, args, 'inversion')
    parent_frame = pd.concat([mother_frame, father_frame])
    # mother_pf = pd.read_csv(args.mpath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # father_pf = pd.read_csv(args.fpath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # raw_pf = pd.concat([mother_pf, father_pf])
    # confident_pf = raw_pf.loc[raw_pf['Confidence'] > 0.3]
    # parent_frame = confident_pf[confident_pf['Type']=='inversion']
    #
    #load reference
    ref_frame = readsmap(args.referencepath, args, 'inversion')
    # raw_rf = pd.read_csv(args.referencepath, sep='\t', comment='#', names=colnames, header=None, skiprows=lambda x: x in [0])
    # confident_rf = raw_rf.loc[raw_rf['Confidence'] > 0.3]
    # del_confident_rf = confident_rf[confident_rf['Type']=='inversion']
    # ref_frame = del_confident_rf

    
    sample_start, sample_end, parent_start, parent_end, ref_start, ref_end = sample_frame.copy(), sample_frame.copy(), parent_frame.copy(), parent_frame.copy(), ref_frame.copy(), ref_frame.copy()

    for df in [sample_start, parent_start, ref_start]: #create an interval for the inversion point
        df['Start'], df['End'], df['Chromosome'] = df.RefStartPos - 100000, df.RefStartPos + 100000, df['RefcontigID1']

    for df in [sample_end, parent_end, ref_end]: #create an interval for the inversion end point
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
    denovo_start = PyRanges(sample_start).overlap(PyRanges(parent_start))
    denovo_end = PyRanges(sample_end).overlap(PyRanges(parent_end))
    denovo_frame = pd.concat([denovo_start.df, denovo_end.df], axis=0)

    if denovo_frame.empty:
        filtered_sample_frame["Found_in_Parent"] = "False"
        calls = filtered_sample_frame.rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()

    else:
        denovo_filtered_sample_frame = pd.merge(filtered_sample_frame, denovo_frame, on=None, how='left', indicator='Found_in_Parent')
        denovo_filtered_sample_frame['Found_in_Parent'] = np.where(denovo_filtered_sample_frame.Found_in_Parent == 'both', True, False)
        calls = denovo_filtered_sample_frame.drop(columns = ['Chromosome', 'Start', 'End']).rename(columns = {'Score':'Phenotype', 'Score2':'Phenotype2'}).drop_duplicates()

    calls.to_csv(args.outputdirectory + '/' + args.sampleID + '_BioNanoInversions.txt', sep='\t', index = False)



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
    parser.add_argument("-g", "--genelist", help="Primary genelist", dest = "genelist", type=str)
    args = parser.parse_args()


    # Actual function
    inversion(args)


if __name__=="__main__":
    main()