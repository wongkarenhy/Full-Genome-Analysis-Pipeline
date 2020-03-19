#!/usr/bin/env python3.6

import pandas as pd
from collections import defaultdict, Counter
import argparse
import sys
import os
import subprocess
import re
import numpy as np
from datetime import datetime
from itertools import chain
from pyranges import PyRanges
from SV_modules import *
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)
pd.options.display.max_rows = 999


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)



def createGeneSyndromeDict(database_df):

    dict = defaultdict(list)
    for var, hpo in database_df.itertuples(index=False):  # var can either be gene or syndrome
        dict[var].append(hpo)

    return(dict)



def createWeightDict(weights):

    try:

        w_df = pd.read_csv(weights, sep = ' ', names=["HPO_id", "weight"], comment = '#')
    except OSError:
        print("Count not open/read the input file:" + weights)
        sys.exit()

    weightDict = dict(zip(w_df.HPO_id, w_df.weight))
    return(weightDict)



def getClinicalPhenome(args):
    # Get the clinical phenome and store as a set
    try:
        clinical_phenome = set(open("./results/" + args.sampleid + "/" + args.sampleid + "_hpo_inexact.txt").read().splitlines())
    except OSError:
        print("Count not open/read the input file:" + "./results/" + args.sampleid + "/" + args.sampleid + "_hpo_inexact.txt")
        sys.exit()

    return(clinical_phenome)



def calculateGeneSumScore(args, hpo_gene_dict, weightDict, clinical_phenome, omim_gene):

    # Go through genes in genelist found in the patients
    try:
        genes = open("./results/" + args.sampleid + "/" + args.sampleid + "_gene_list.txt", 'r')
    except OSError:
        print("Count not open/read the input file:" + "./results/" + args.sampleid + "/" + args.sampleid + "_gene_list.txt")
        sys.exit()

    with genes:

        gene = genes.read().splitlines()
        gene_sum_score = 0
        gene_score_result = pd.DataFrame(columns=['gene', 'score'])

        for query in gene:
            #print(query)
            hpo_pheno = set(hpo_gene_dict[query]) # To get the phenotypic features for a given gene
            overlap = hpo_pheno.intersection(clinical_phenome) # overlap all the phenotypic features with the clinical phenomes

            for term in overlap:
                gene_sum_score += weightDict[term]

            gene_score_result = gene_score_result.append({'gene':query, 'score':gene_sum_score}, ignore_index=True)

        gene_score_result_r = gene_score_result.iloc[::-1]
        gene_score_result_r = pd.concat([gene_score_result_r, omim_gene])
        gene_score_result_r = normalizeRawScore(args, gene_score_result_r, 'gene')

    return(gene_score_result_r)



def getParentsGeno(filtered_intervar, inheritance_mode, ov_allele):

    # Create two new columns and initialize to 0
    filtered_intervar[inheritance_mode] = 0
    filtered_intervar = filtered_intervar.reset_index(drop=True)

    for idx, row in enumerate(filtered_intervar.itertuples(index=False)):
        if int(getattr(row, 'Start')) in set(ov_allele['Start']):
            parents_geno = ov_allele.loc[ov_allele['Start']==getattr(row,'Start'),'geno'].item()
            filtered_intervar.loc[idx, inheritance_mode] = parents_geno

    return(filtered_intervar)




def rerankSmallVariant(df):

    df['Clinvar_idx'] = df.Clinvar.str[9:-1]
    df['InterVar_idx'] = df.InterVar_InterVarandEvidence.str[10:].str.split('PVS1').str[0]
    df[['Clinvar_idx', 'InterVar_idx']] = df[['Clinvar_idx', 'InterVar_idx']].apply(lambda x:x.astype(str).str.lower())
    df['Clinvar_score'], df['InterVar_score'] = 3, 3

    # Calculate Clinvar score
    df.loc[(df['Clinvar_idx'].str.contains('benign')), 'Clinvar_score'] = 1
    df.loc[((df['Clinvar_idx'].str.contains('benign')) & (df['Clinvar_idx'].str.contains('likely'))), 'Clinvar_score'] = 2
    df.loc[(df['Clinvar_idx'].str.contains('pathogenic')), 'Clinvar_score'] = 5
    df.loc[((df['Clinvar_idx'].str.contains('pathogenic')) & (df['Clinvar_idx'].str.contains('likely'))), 'Clinvar_score'] = 4
    df.loc[(df['Clinvar_idx'].str.contains('conflicting')), 'Clinvar_score'] = 3

    # Calculate Intervar score
    df.loc[(df['InterVar_idx'].str.contains('benign')), 'InterVar_score'] = 1
    df.loc[((df['InterVar_idx'].str.contains('benign')) & (df['InterVar_idx'].str.contains('likely'))), 'InterVar_score'] = 2
    df.loc[(df['InterVar_idx'].str.contains('pathogenic')), 'InterVar_score'] = 5
    df.loc[((df['InterVar_idx'].str.contains('pathogenic')) & (df['InterVar_idx'].str.contains('likely'))), 'InterVar_score'] = 4
    # Add them up
    df['Patho_score'] = df['Clinvar_score'] + df['InterVar_score']

    # Sort by the total patho_score
    df = df.sort_values(by=['Patho_score', 'score'], ascending=False)
    df = df.drop(['Clinvar_idx', 'InterVar_idx', 'Clinvar_score', 'InterVar_score', 'Patho_score'], axis=1)

    return df



def smallVariantGeneOverlapCheckInheritance(args, smallVariantFile, interVarFinalFile, gene_score_result_r, famid):

    # Overlap gene_score_result_r with small variants genes found in the proband
    gene_score_result_r = gene_score_result_r[gene_score_result_r.gene.isin(smallVariantFile.gene)]


    # Subset the intervar files further to store entries relevant to these set of genes
    filtered_intervar = pd.merge(interVarFinalFile, gene_score_result_r, left_on='Ref_Gene', right_on='gene',how='inner')

    # Remove common artifacts
    try:
        artifacts = pd.read_csv("./common_artifacts_20.txt", names = ["gene"])
    except OSError:
        print("Could not open/read the input file: common_artifacts_20.txt")
        sys.exit()

    filtered_intervar = filtered_intervar.loc[~filtered_intervar['Ref_Gene'].isin(artifacts['gene'])]

    # Create a bed file and write it out
    pd.DataFrame(filtered_intervar).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_smallVariant_candidates.txt', index=False, sep='\t',header=False)  # Write out a subset of the variant first
    filtered_intervar_bed = filtered_intervar[['Chr', 'Start', 'End']]
    filtered_intervar_bed.loc[:,'Chr'] = 'chr' + filtered_intervar_bed.loc[:,'Chr'].astype(str)
    filtered_intervar_bed.loc[:,'Start'] -= 1
    pd.DataFrame(filtered_intervar_bed).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_target.bed', index=False, sep='\t', header=False)

    if not args.singleton:

        # Get overlapping variants from the parents so we know which variants are inherited
        print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Comparing small variants (SNPs/indels) inheritance')
        cmd1 = "bcftools view -R ./results/" + args.sampleid + "/" + args.sampleid + "_target.bed /media/KwokRaid04/CIAPM/CIAPM_longranger/BC0" + famid + "01_longranger/outs/phased_variants.vcf.gz > ./results/" + args.sampleid + "/" + args.sampleid + "_paternal_inherited_smallVariants.vcf"
        cmd2 = "bcftools view -R ./results/" + args.sampleid + "/" + args.sampleid + "_target.bed /media/KwokRaid04/CIAPM/CIAPM_longranger/BC0" + famid + "02_longranger/outs/phased_variants.vcf.gz > ./results/" + args.sampleid + "/" + args.sampleid + "_maternal_inherited_smallVariants.vcf"
        cmds = [cmd1, cmd2]
        for cmd in cmds:
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            if p.returncode != 0:
                raise Exception(stderr)

        # Go through every row in filtered_intervar and see if the same variant is found in either of the parents
        # We will only compare allele start position (we always assume the alt allele is the same)
        try:
            paternal_ov_allele = pd.read_csv("./results/" + args.sampleid + "/" + args.sampleid + "_paternal_inherited_smallVariants.vcf", sep='\t',usecols=[1,9], names=["Start", "geno"], comment='#')
            paternal_ov_allele['geno'] = paternal_ov_allele['geno'].str[:1].astype(int) + paternal_ov_allele['geno'].str[2:3].astype(int)
        except OSError:
            print("Could not open/read the input file: ./results/" + args.sampleid + "/" + args.sampleid + "_paternal_inherited_smallVariants.vcf")
            sys.exit()
        try:
            maternal_ov_allele = pd.read_csv("./results/" + args.sampleid + "/" + args.sampleid + "_maternal_inherited_smallVariants.vcf", sep='\t',usecols=[1,9], names=["Start", "geno"], comment='#')
            maternal_ov_allele['geno'] = maternal_ov_allele['geno'].str[:1].astype(int) + maternal_ov_allele['geno'].str[2:3].astype(int)
        except OSError:
            print(
                "Could not open/read the input file: ./results/" + args.sampleid + "/" + args.sampleid + "_maternal_inherited_smallVariants.vcf")
            sys.exit()


        filtered_intervar = getParentsGeno(filtered_intervar, 'paternal', paternal_ov_allele)
        filtered_intervar = getParentsGeno(filtered_intervar, 'maternal', maternal_ov_allele)

        # Rerank variants based on reported or predicted pathogeneicity
        filtered_intervar = rerankSmallVariant(filtered_intervar)

        # Divide the dataset into recessive, dominant, de novo, compound het
        ## Recessive
        recessive = filtered_intervar[(filtered_intervar['paternal'] == 1) & (filtered_intervar['maternal'] == 1) & (filtered_intervar['Otherinfo'] == 'hom')]
        ## Dominant
        dominant_inherited = filtered_intervar[((filtered_intervar['paternal'] == 1) & (filtered_intervar['maternal'] == 0)) | ((filtered_intervar['maternal'] == 1) & (filtered_intervar['paternal'] == 0))]
        ## De novo
        denovo = filtered_intervar[(filtered_intervar['paternal'] == 0) & (filtered_intervar['maternal'] == 0)]
        #Compound het
        filtered_intervar_compoundhet = filtered_intervar[(filtered_intervar['Otherinfo'] == 'het')]
        filtered_intervar_compoundhet = filtered_intervar_compoundhet[(filtered_intervar_compoundhet['maternal'] != 2) & (filtered_intervar_compoundhet['paternal'] != 2) & ((filtered_intervar_compoundhet['paternal'] == 1) & (filtered_intervar_compoundhet['maternal'] == 0)) | ((filtered_intervar_compoundhet['maternal'] == 1) & (filtered_intervar_compoundhet['paternal'] == 0)) | ((filtered_intervar_compoundhet['maternal'] == 0) & (filtered_intervar_compoundhet['paternal'] == 0))]
        count = Counter(filtered_intervar_compoundhet['Ref_Gene'])
        compoundhet_genes = [x for x, cnt in count.items() if cnt > 1]
        compoundhet = filtered_intervar_compoundhet[filtered_intervar_compoundhet['Ref_Gene'].isin(compoundhet_genes)]

        discard = []
        for gene in compoundhet_genes:

            df = compoundhet[compoundhet['Ref_Gene'].str.contains(gene)]
            row_count = len(df.index)
            col_list = ['paternal', 'maternal']
            res = df[col_list].sum(axis=0)
            if ((res[0] == 0) & (res[1] == row_count)) or (res[1] == 0 & (res[0] == row_count)):
                discard.append(gene)

        compoundhet = compoundhet[~compoundhet['Ref_Gene'].isin(discard)]


        # Print all the variants according to inheritance mode
        # Recessive
        pd.DataFrame(recessive).to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_smallVariants_recessive_candidates.txt', index=False, sep='\t', header=True)
        # Dominant
        pd.DataFrame(dominant_inherited).to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_dominant_inherited_smallVariants_candidates.txt', index=False, sep='\t', header=True)
        # De novo
        pd.DataFrame(denovo).to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_smallVariants_denovo_candidates.txt', index=False, sep='\t', header=True)
        # Compound het
        pd.DataFrame(compoundhet).to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_smallVariants_compoundhet_candidates.txt', index=False, sep='\t', header=True)

    # All
    filtered_intervar = rerankSmallVariant(filtered_intervar)
    pd.DataFrame(filtered_intervar).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_smallVariants_ALL_candidates.txt', index=False, sep='\t', header=True)

    if not args.singleton:
        # We want to return everything except recessive variants
        filtered_intervar = filtered_intervar.loc[~filtered_intervar['Start'].isin(recessive['Start'])] # don't have recessive if singleton

    return filtered_intervar



def differentialDiangosis(hpo_syndrome_dict, weightSyndromeDict, clinical_phenome, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup,hpo_syndromes_mim_df):

    syndrome_score_result = pd.DataFrame(columns=['syndrome', 'score'])


    # Check every syndrome and its overlapping hpo terms
    for syndrome in hpo_syndrome_dict:

        hpo_terms = set(hpo_syndrome_dict[syndrome])

        score = 0
        for term in hpo_terms:

            if term in clinical_phenome:
                score += weightSyndromeDict[term]

        if score != 0:
            syndrome_score_result = syndrome_score_result.append({'syndrome': syndrome, 'score': score}, ignore_index=True)

    syndrome_score_result_r = syndrome_score_result.sort_values(by='score', ascending=False)
    syndrome_score_result_r['syndrome'] = syndrome_score_result_r['syndrome'].str.upper()

    # Add a normalized score column
    syndrome_score_result_r = normalizeRawScore(args, syndrome_score_result_r, 'syndrome')

    # Specifically look for deletion/duplication syndrome
    delDupSyndrome(syndrome_score_result_r, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup, hpo_syndromes_mim_df)

    return(syndrome_score_result_r)




def findGenomicLocation(cytoband_key, cytobandDict):
    print(cytoband_key)
    keys = [key for key in cytobandDict if key.startswith(cytoband_key)]
    print(keys)
    if len(keys)==0:
        cytoband_key = cytoband_key[:-1]
        keys = [key for key in cytobandDict if key.startswith(cytoband_key)]

    genomic_coords_list = []

    for key in keys:
        genomic_coords_list.append(str(cytobandDict[key]).split('-'))
    print(genomic_coords_list)
    genomic_coords_list = list(chain.from_iterable(genomic_coords_list))
    min_coords = min(genomic_coords_list)
    max_coords = max(genomic_coords_list)

    genomic_range = str(min_coords) + '-' + str(max_coords)

    return genomic_range





def parseSyndromeNameToCytoband(df, cytobandDict, type, hpo_syndromes_mim_df,args):

    if type=='deldup':
        df['cytoband'] = float('Nan')
        regex = r'((^|\W)[0-9XY]{1,2}[PQ]{1}[\w\\.\\-]{1,15}[\s$])'

        for index, row in df.iterrows():
            m = re.search(regex, str(row))
            if m is not None:
                df.loc[index, 'cytoband'] = m.group(1)

        df.dropna(subset=['cytoband'], inplace=True)
        if df.empty:  # df can be empty after dropping NA
            return pd.DataFrame()

    if type=='all':
        df = df.merge(hpo_syndromes_mim_df, on=['syndrome'])

        try:
            morbid = pd.read_csv(args.workdir + '/morbidmap.txt', sep='\t', usecols=[2, 3], names=["MIM", "cytoband"], comment='#')
            df = df.merge(morbid, on='MIM')
            df = df.loc[~df['cytoband'].astype(str).str.contains("Chr")]
            end_string = ('p','q')
            df = df.loc[~df['cytoband'].str.endswith(end_string)] #Remove cytoband entries that span the whole chromosomal arm like 2p
            print(df)

        except OSError:
            print("Could not open/read the input file: " + args.workdir + '/morbidmap.txt')
            sys.exit()

    df['cytoband'] = df['cytoband'].astype(str).str.lower()
    df['cytoband'] = df['cytoband'].str.replace('x', 'X')
    df['cytoband'] = df['cytoband'].str.replace('y', 'Y')
    df['cytoband'] = df['cytoband'].str.strip('\(\)')
    df[['Chromosome', 'discard']] = df.cytoband.str.split('p|q', 1, expand=True)
    df = df.drop('discard', axis=1)
    if df.cytoband.str.contains('-').any():
        df[['cytoband_start', 'cytoband_stop']] = df.cytoband.str.split('-', expand=True)
    else:
        df['cytoband_start'] = df.cytoband
        df['cytoband_stop'] = None
    df['arm'] = np.where(df['cytoband_start'].str.contains('p'), 'p', 'q')

    df['cytoband_stop'] = np.where(df['cytoband_start'].str.count('p|q')>1, df['arm'] + df['cytoband_start'].str.split('p|q').str[2], df['cytoband_stop'])
    df['cytoband_start'] = np.where(df['cytoband_start'].str.count('p|q')>1, df['cytoband_start'].str.split('p|q').str[0] + df['arm'] + df['cytoband_start'].str.split('p|q').str[1], df['cytoband_start'])


    for idx, row in df.iterrows():
        cytoband_start_key = row['cytoband_start'].replace(" ","")
        if cytoband_start_key in cytobandDict:
            coords_start = cytobandDict[cytoband_start_key]

        else:
            genomic_range = findGenomicLocation(cytoband_start_key, cytobandDict)
            coords_start = genomic_range


        if row['cytoband_stop'] is not None: # Fix cytoband_stop column for quick cytobandDict lookup
            current_chr = np.where(('p' in str(row['cytoband_stop'])) or ('q' in str(row['cytoband_stop'])), str(row['Chromosome']), str(row['Chromosome']) + str(row['arm']))
            edited_cytoband_stop = str(current_chr) + row['cytoband_stop']
            edited_cytoband_stop = edited_cytoband_stop.replace(" ", "")
            df.at[idx, 'cytoband_stop'] = edited_cytoband_stop

            if edited_cytoband_stop in cytobandDict:
                coords_stop = cytobandDict[edited_cytoband_stop]
            else:
                genomic_range = findGenomicLocation(edited_cytoband_stop, cytobandDict)
                coords_stop = genomic_range

            # New coords will be the the beginning of coords_start and end of coords_stop
            df.at[idx, 'Start'] = coords_start.split('-')[0]
            df.at[idx, 'End'] = coords_stop.split('-')[1]

        else:
            df.at[idx, 'Start'] = coords_start.split('-')[0]
            df.at[idx, 'End'] = coords_start.split('-')[1]

    return df




def createCytobandDict(args):

    try:
        cyto = pd.read_csv(args.workdir + '/cytoband.txt', sep = '\t', names=["cytoband", "coords"], comment = '#')
    except OSError:
        print("Count not open/read the input file:" + args.workdir + '/cytoband.txt')
        sys.exit()

    cytobandDict = dict(zip(cyto.cytoband, cyto.coords))

    return(cytobandDict)




def delDupSyndrome(syndrome_score_result_r, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup, hpo_syndromes_mim_df):

    #print(syndrome_score_result_r)
    syndrome_score_result_r.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_syndrome_score_result_r.txt', sep='\t', index=False)
    # Create cytoband <-> genomic coordinates dict
    cytobandDict = createCytobandDict(args)

    del_cond = syndrome_score_result_r['syndrome'].str.contains('DELETION')
    dup_cond = syndrome_score_result_r['syndrome'].str.contains('DUPLICATION')

    del_df = syndrome_score_result_r[del_cond]
    dup_df = syndrome_score_result_r[dup_cond]

    del_df = parseSyndromeNameToCytoband(del_df, cytobandDict,'deldup',hpo_syndromes_mim_df, args)
    dup_df = parseSyndromeNameToCytoband(dup_df, cytobandDict,'deldup',hpo_syndromes_mim_df, args)

    all_omim_syndromes = parseSyndromeNameToCytoband(syndrome_score_result_r, cytobandDict,'all', hpo_syndromes_mim_df, args)

    # print(dup_df)
    #
    # if del_df.empty:
    #     del_df = None
    # else:
    #     del_df = parseSyndromeNameToCytoband(del_df, cytobandDict)
    # if dup_df.empty:
    #     dup_df = None
    # else:
    #     dup_df = parseSyndromeNameToCytoband(dup_df, cytobandDict)

    if args.bionano:
        if not args.singleton:
            cols = ['Chromosome', 'Start', 'End', 'SmapEntryID', 'Confidence', 'Type', 'Zygosity', 'Genotype', 'SV_size', 'Found_in_Father', 'Found_in_Mother', 'syndrome', 'cytoband', 'score', 'normalized_score']
        else:
            cols = ['Chromosome', 'Start', 'End', 'SmapEntryID', 'Confidence', 'Type', 'Zygosity', 'Genotype', 'SV_size', 'syndrome', 'cytoband', 'score', 'normalized_score']

        # Overlap with del/dup syndromes
        if cyto_BN_dup is not None: # It can be None because old Bionano pipeline doesn't call duplications...
            # dup_df.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_input1.txt', sep='\t', index=False)
            # cyto_BN_dup.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_input2.txt', sep='\t',index=False)
            # cyto_10x_dup_largeSV.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_input2.txt', sep='\t',index=False)

            overlap_dup_BN = delDupSyndromeSVOverlap(dup_df, cyto_BN_dup, cols)
            overlap_dup_BN.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_duplication_syndrome.txt', sep='\t', index=False)
        else:
            overlap_dup_BN = None
            pd.DataFrame().to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_duplication_syndrome.txt', sep='\t', index=False)

        overlap_del_BN = delDupSyndromeSVOverlap(del_df, cyto_BN_del, cols)
        overlap_del_BN.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_deletion_syndrome.txt', sep='\t', index=False)

        all_BN = pd.concat([cyto_BN_dup, cyto_BN_del])
        overlap_all_BN = delDupSyndromeSVOverlap(all_omim_syndromes, all_BN, cols)
        overlap_all_BN.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_all_syndrome.txt', sep='\t', index=False)

    if args.linkedreadSV:

        if not args.singleton:
            cols = ['Chromosome', 'Start', 'End', 'ID', 'REF', 'ALT_1', 'QUAL', 'FILTER_PASS', 'SVLEN', 'Found_in_Father', 'Found_in_Mother', 'syndrome', 'cytoband', 'score', 'normalized_score']
        else:
            cols = ['Chromosome', 'Start', 'End', 'ID', 'REF', 'ALT_1', 'QUAL', 'FILTER_PASS', 'SVLEN', 'syndrome', 'cytoband', 'score', 'normalized_score']

        # dup_df.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_input1.txt', sep='\t', index=False)
        # cyto_10x_dup_largeSV.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_input2.txt', sep='\t', index=False)

        overlap_dup_largeSV_10x = delDupSyndromeSVOverlap(dup_df, cyto_10x_dup_largeSV, cols)
        overlap_dup_largeSV_10x.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_duplication_largeSV_syndrome.txt', sep='\t', index=False)

        overlap_del_largeSV_10x = delDupSyndromeSVOverlap(del_df, cyto_10x_del_largeSV, cols)
        overlap_del_largeSV_10x.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_deletion_largeSV_syndrome.txt', sep='\t', index=False)


        overlap_del_10x = delDupSyndromeSVOverlap(del_df, cyto_10x_del, cols)
        overlap_del_10x.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_deletion_syndrome.txt', sep='\t', index=False)

        all_10x = pd.concat([cyto_10x_dup_largeSV, cyto_10x_del_largeSV, cyto_10x_del], ignore_index=True)
        overlap_all_10x = delDupSyndromeSVOverlap(all_omim_syndromes, all_10x, cols)
        overlap_all_10x.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_all_syndrome.txt', sep='\t', index=False)


    if args.linkedreadSV and args.bionano:

        if not args.singleton:
            cols = ['Chromosome', 'Start', 'End', 'ID', 'REF', 'ALT_1', 'QUAL', 'FILTER_PASS', 'SVLEN', 'Found_in_Father', 'Found_in_Mother', 'syndrome', 'cytoband', 'SmapEntryID', 'Confidence', 'Type', 'Zygosity', 'Genotype', 'SV_size', 'Found_in_Father_b', 'Found_in_Mother_b', 'score', 'normalized_score',]
        else:
            cols = ['Chromosome', 'Start', 'End', 'ID', 'REF', 'ALT_1', 'QUAL', 'FILTER_PASS', 'SVLEN', 'syndrome', 'cytoband', 'SmapEntryID', 'Confidence', 'Type', 'Zygosity', 'Genotype', 'SV_size', 'score', 'normalized_score',]

        # syndrome appearing in both 10x and bionano --> confident set
        ## for duplications
        if ((overlap_dup_BN is not None) and (not overlap_dup_BN.empty) and (not overlap_dup_largeSV_10x.empty)):
            overlap_dup_largeSV_10x = overlap_dup_largeSV_10x.loc[overlap_dup_largeSV_10x['SVLEN'] >= 1000]
            confident_dup_syndrome = delDupSyndromeSVOverlap(overlap_dup_largeSV_10x, overlap_dup_BN, cols)

            if not confident_dup_syndrome.empty:
                confident_dup_syndrome.to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_duplication_syndrome.txt', sep='\t',index=False)
            else: # Write an empty dataframe
                pd.DataFrame().to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_duplication_syndrome.txt', sep='\t',index=False)

        ## for deletions
        del_10x = pd.concat([overlap_del_largeSV_10x, overlap_del_10x])
        if ((not overlap_del_BN.empty) and (not del_10x.empty)):
            del_10x = del_10x.loc[del_10x['SVLEN'] <= (-1000)]
            confidnet_del_syndrome = delDupSyndromeSVOverlap(del_10x, overlap_del_BN, cols)

            #confidnet_del_syndrome = pd.merge(del_10x, overlap_del_BN, on='syndrome', how='inner')
            if not confidnet_del_syndrome.empty:
                confidnet_del_syndrome.to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_deletion_syndrome.txt', sep='\t',index=False)
            else:
                pd.DataFrame().to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_deletion_syndrome.txt', sep='\t',index=False)

        # for all omim syndromes
        if ((not overlap_all_BN.empty) and (not overlap_all_10x.empty)):
            overlap_all_10x = overlap_all_10x.loc[(overlap_all_10x['SVLEN'] <= (-1000)) | (overlap_all_10x['SVLEN'] >=1000)]
            confident_all_syndrome = delDupSyndromeSVOverlap(overlap_all_10x, overlap_all_BN, cols)

            if not confident_all_syndrome.empty:
                confident_all_syndrome.to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_all_syndrome.txt', sep='\t',index=False)
            else:
                pd.DataFrame().to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_all_syndrome.txt', sep='\t',index=False)




def delDupSyndromeSVOverlap(del_df, cyto_BN_del, cols):

    if del_df.empty:
        return pd.DataFrame()

    del_df['Chromosome'] = del_df['Chromosome'].str.strip()
    del_df.dropna( inplace=True)

    overlap_del_BN = PyRanges(cyto_BN_del).join(PyRanges(del_df))

    if not overlap_del_BN.df.empty:
        overlap_del_BN = overlap_del_BN.df
        overlap_del_BN['overlap_len'] = np.maximum(0, np.minimum(overlap_del_BN.End, overlap_del_BN.End_b) - np.maximum(overlap_del_BN.Start,overlap_del_BN.Start_b))
        #overlap_del_BN = overlap_del_BN.drop(like="_b")
        overlap_del_BN = overlap_del_BN.sort_values(by='score', ascending=False)
        overlap_del_BN = overlap_del_BN.loc[overlap_del_BN['overlap_len'] > 0]
        # print(overlap_del_BN)
        #overlap_del_BN = overlap_del_BN.df.sort_values(by='score', ascending=False)
        # Rearrange the column
        overlap_del_BN = overlap_del_BN[cols].drop_duplicates()
        return overlap_del_BN
    else:
        return overlap_del_BN.df





def normalizeRawScore(args, raw_score, mode):

    # Normalize all the scores to 1-100
    max_score = max(raw_score['score'])
    raw_score.loc[:,'normalized_score'] = raw_score.loc[:,'score']/max_score * 100

    return(raw_score)



def compileControlFiles(control_files_path, famid):

    full_paths = []
    for path in control_files_path:
        control_files = os.listdir(path)
        for file in control_files:
            if not (re.match('BC0..0[34]{1}', file) or re.match(rf"BC0{famid}..", file)):  # Discard trio of interest and all probands
                full_paths.append(os.path.join(path, file))
                full_paths.append(os.path.join(path, file))

    return full_paths





def bionanoSV(args, famid, gene_score_result_r, all_small_variants):

    # # Generate controls files (1KGP BN samples + CIAPM parents (excluding  parents of the proband of interest)
    # print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Generating bionano control file...')
    # control_files_path = [args.workdir + "/bionano_sv/controls/DLE", args.workdir + "/bionano_sv/controls/BspQI", args.workdir + "/bionano_sv/cases/DLE", args.workdir + "/bionano_sv/cases/BspQI"]
    # full_paths = compileControlFiles(control_files_path, famid)
    #
    # ## Write an empty file
    # with open(args.workdir + "/results/" + args.sampleid + "/bionano_control.smap.gz", 'w'):  # So it will overwrite the old file
    #     pass
    #
    # for path in full_paths:
    #     cmd = "cat " + path + "/exp_refineFinal1_merged_filter.smap | gzip >> " + args.workdir + "/results/" + args.sampleid + "/bionano_control.smap.gz"
    #     p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #     stdout, stderr = p.communicate()
    #     if p.returncode != 0:
    #         raise Exception(stderr)

    # Create a BN arg object
    BN_args = Namespace(sampleID = args.sampleid,
                        samplepath = args.workdir + "/bionano_sv/cases/" + args.enzyme + "/" + args.sampleid + "/exp_refineFinal1_merged_filter.smap",
                        fpath = args.workdir + "/bionano_sv/cases/" + args.enzyme + "/BC0" + famid + "01/exp_refineFinal1_merged_filter.smap",
                        mpath = args.workdir + "/bionano_sv/cases/" + args.enzyme + "/BC0" + famid + "02/exp_refineFinal1_merged_filter.smap",
                        referencepath = args.workdir + "/results/" + args.sampleid + "/bionano_control.smap.gz",
                        outputdirectory = args.workdir + '/results/' + args.sampleid,
                        exons = args.workdir + '/annotatedExon.bed',
                        genes=args.workdir + '/annotatedGene.bed',
                        genelist = gene_score_result_r,
                        singleton = args.singleton)


    # # Call bionano translocation
    # print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano translocations on ' + args.sampleid + '...')
    # BN_translocation(BN_args)

    # Call bionano deletion
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano deletions on ' + args.sampleid + '...')
    cyto_BN_del, exon_calls_BN_del = BN_deletion(BN_args)

    # Call bionano insertion
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano insertions on ' + args.sampleid + '...')
    BN_insertion(BN_args)

    # Call bionano duplications
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano duplications on ' + args.sampleid + '...')
    cyto_BN_dup, exon_calls_BN_dup = BN_duplication(BN_args)

    # # Call bionano inversions
    # print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano inversions on ' + args.sampleid + '...')
    # BN_inversion(BN_args)

    # Check potential compoundhets with SNPs and indels
    BN_exons = pd.concat([exon_calls_BN_del, exon_calls_BN_dup])
    if BN_exons.empty:
        pd.DataFrame().to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_Bionano_SV_SNPsIndels_compoundhet_candidates.txt', sep='\t', index=False)
    else:
        BN_exons = pd.merge(BN_exons, all_small_variants, left_on='gene', right_on='Ref_Gene', how='inner')
        BN_exons.to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_Bionano_SV_SNPsIndels_compoundhet_candidates.txt', sep='\t', index=False)

    return cyto_BN_del, cyto_BN_dup, exon_calls_BN_del, exon_calls_BN_dup




def linkedreadSV(args, famid, gene_score_result_r, all_small_variants):

    # Need to generate a reference file for all the medium size deletions
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Generating linked-reads control files...')
    control_files_path = [args.workdir + "/linkedRead_sv/controls", args.workdir + "/linkedRead_sv/cases"]
    full_paths = compileControlFiles(control_files_path, famid)
    ## Write an empty file
    with open(args.workdir + "/results/" + args.sampleid + "/10x_del_control.vcf.gz",'w'):  # So it will overwrite the old file
        pass

    for path in full_paths:
        cmd = "zcat " + path + "/dels.vcf.gz | gzip >> " + args.workdir + "/results/" + args.sampleid + "/10x_del_control.vcf.gz"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise Exception(stderr)


    # Need to generate another reference file for large SVs
    with open(args.workdir + "/results/" + args.sampleid + "/10x_largeSV_control.vcf.gz",'w'):  # So it will overwrite the old file
        pass

    for path in full_paths:
        cmd = "zcat " + path + "/large_svs.vcf.gz | gzip >> " + args.workdir + "/results/" + args.sampleid + "/10x_largeSV_control.vcf.gz"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise Exception(stderr)

    tenx_args_del = Namespace(sampleID = args.sampleid,
                        samplepath = args.workdir + "/linkedRead_sv/cases/" + args.sampleid + "/dels.vcf.gz",
                        fpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "01/dels.vcf.gz",
                        mpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "02/dels.vcf.gz",
                        referencepath = args.workdir + "/results/" + args.sampleid + "/10x_del_control.vcf.gz",
                        outputdirectory = args.workdir + '/results/' + args.sampleid,
                        exons = args.workdir + '/annotatedExon.bed',
                        genes = args.workdir + '/annotatedGene.bed',
                        genelist = gene_score_result_r,
                        singleton = args.singleton)

    tenx_args_largeSV = Namespace(sampleID = args.sampleid,
                        samplepath = args.workdir + "/linkedRead_sv/cases/" + args.sampleid + "/large_svs.vcf.gz",
                        fpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "01/large_svs.vcf.gz",
                        mpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "02/large_svs.vcf.gz",
                        referencepath = args.workdir + "/results/" + args.sampleid + "/10x_largeSV_control.vcf.gz",
                        outputdirectory = args.workdir + '/results/' + args.sampleid,
                        exons = args.workdir + '/annotatedExon.bed',
                        genes=args.workdir + '/annotatedGene.bed',
                        genelist = gene_score_result_r,
                        singleton = args.singleton)

    # Call medium size deletions
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads medium deletions on ' + args.sampleid + '...')
    cyto_10x_del, exon_calls_10x_del = tenxdeletions(tenx_args_del)

    # Call large deletions
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large deletions on ' + args.sampleid + '...')
    cyto_10x_del_largeSV, exon_calls_10x_largeSV_del = tenxlargesvdeletions(tenx_args_largeSV)

    # Call large duplications
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large duplications on ' + args.sampleid + '...')
    cyto_10x_dup_largeSV, exon_calls_10x_largeSV_dup = tenxlargesvduplications(tenx_args_largeSV)

    # # Call large inversions
    # print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large inversions on ' + args.sampleid + '...')
    # tenxlargesvinversions(tenx_args_largeSV)
    #
    # # Call large breakends
    # print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large breakends on ' + args.sampleid + '...')
    # tenxlargesvbreakends(tenx_args_largeSV)
    #
    # # Call large unknwon calls
    # print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large unknown on ' + args.sampleid + '...')
    # tenxlargesvunknown(tenx_args_largeSV)

    # Check potential compoundhets with SNPs and indels
    tenx_exons = pd.concat([exon_calls_10x_del, exon_calls_10x_largeSV_del, exon_calls_10x_largeSV_dup])
    if tenx_exons.empty:
        pd.DataFrame().to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_10x_SV_SNPsIndels_compoundhet_candidates.txt', sep='\t', index=False)
    else:
        tenx_exons = pd.merge(tenx_exons, all_small_variants, left_on='gene', right_on='Ref_Gene', how='inner')
        tenx_exons.to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_10x_SV_SNPsIndels_compoundhet_candidates.txt', sep='\t', index=False)

    return cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, exon_calls_10x_del, exon_calls_10x_largeSV_del, exon_calls_10x_largeSV_dup



def pyrangeJoin(df1_10x, df2_BN):

    if df1_10x.empty or df2_BN.empty:
        return pd.DataFrame()

    df1_10x['Chromosome'], df1_10x['Start'], df1_10x['End'] = [df1_10x['CHROM'], df1_10x['POS'], df1_10x['END']]
    df2_BN['Chromosome'], df2_BN['Start'], df2_BN['End'] = [df2_BN['RefcontigID1'], df2_BN['RefStartPos'], df2_BN['RefEndPos']]
    overlap = PyRanges(df1_10x).join(PyRanges(df2_BN))

    #print(overlap)
    if not overlap.df.empty:
        overlap = overlap.df
        overlap['overlap_len'] = np.maximum(0, np.minimum(overlap.End, overlap.End_b) - np.maximum(overlap.Start,overlap.Start_b))
        #overlap = overlap.drop(like="_b")
        overlap = overlap.drop(['Chromosome', 'Start', 'End'], axis = 1)
        overlap = overlap.loc[overlap['overlap_len'] > 0]
        return overlap
    else:
        return overlap.df



def findConfDelDup(args, exon_calls_10x_del, exon_calls_10x_largeSV_del, exon_calls_10x_largeSV_dup, exon_calls_BN_del, exon_calls_BN_dup):

    tenx_del = pd.concat([exon_calls_10x_del, exon_calls_10x_largeSV_del])
    overlap_del_10x_BN = pyrangeJoin(tenx_del, exon_calls_BN_del)
    overlap_del_10x_BN.to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_deletion_exons.txt', sep='\t', index=False)

    if exon_calls_BN_dup is not None:  # some bionano assemblies were generated with old pipelines
        overlap_dup_10x_BN = pyrangeJoin(exon_calls_10x_largeSV_dup, exon_calls_BN_dup)
        overlap_dup_10x_BN.to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_duplication_exons.txt', sep='\t', index=False)
    else:
        pd.DataFrame().to_csv('./results/' + args.sampleid + "/confident_set/" + args.sampleid + '_confident_duplication_exons.txt', sep='\t', index=False)




def main():

    # Parse argument
    parser = argparse.ArgumentParser(description="This software ranks genes based on the clinical phenome.")
    parser.add_argument("-s", "--sampleid",help="Sample ID",dest="sampleid", type=str, required = True)
    parser.add_argument("-w", "--workdir", help="This is the base work directory.", dest="workdir", type=str, required = True)
    parser.add_argument("-d", "--database", help="Path to HPO database", dest="database", type=str, required = True)
    parser.add_argument("-i", "--intervar", help="Path to InterVar output folder", dest="intervar", type=str, required = True)
    parser.add_argument("-b", "--bionano", help="Set this flag to evaluate bionano SVs.", dest="bionano", action='store_true')
    parser.add_argument("-l", "--linkedreadSV", help="Set this flag to evaluate linkedread SVs.", dest="linkedreadSV", action='store_true')
    parser.add_argument("-e", "--enzyme", help="Bionano enzyme used (BspQI or DLE). Only set this flag if -b is set", dest="enzyme", type=str)
    parser.add_argument("-S", help="Set this flag if this is a singleton case", dest="singleton", action='store_true')
    args = parser.parse_args()

    # Change work dir
    os.chdir(args.workdir)

    # Define variables
    ## Read the database files
    hpo_genes = args.database + "/genes_to_phenotype.txt"
    hpo_syndromes = args.database + "/phenotype_annotation.tab"
    omim_gene = args.workdir + "/annotatedGene.bed"
    smallVariantFileName = args.intervar + "/example/" + args.sampleid + "_smallVariant_geneList.txt"
    interVarFinalFileName = args.intervar + "/example/" + args.sampleid +  ".hg38_multianno.txt.intervar.FINAL"

    try:
        hpo_genes_df = pd.read_csv(hpo_genes, sep='\t', usecols=[1, 2], names=["gene_name", "HPO_id"], comment='#')
    except OSError:
        print("Could not open/read the input file: " + hpo_genes)
        sys.exit()

    try:
        hpo_syndromes_df = pd.read_csv(hpo_syndromes, sep='\t', usecols=[2, 4], names=["syndrome", "HPO_id"], comment='#')
        hpo_syndromes_mim_df = pd.read_csv(hpo_syndromes, sep='\t', usecols=[1, 2], names=["MIM","syndrome"],comment='#')
        hpo_syndromes_mim_df = hpo_syndromes_mim_df.drop_duplicates()
        hpo_syndromes_mim_df['syndrome'] = hpo_syndromes_mim_df['syndrome'].str.upper()
    except OSError:
        print("Could not open/read the input file: " + hpo_syndromes)
        sys.exit()

    try:
        smallVariantFile = pd.read_csv(smallVariantFileName, names=["gene"])
    except OSError:
        print ("Could not open/read the input file: " + smallVariantFileName)

    try:
        interVarFinalFile = pd.read_csv(interVarFinalFileName, sep='\t', index_col = False, names=["Chr", "Start", "End", "Ref", "Alt", "Ref_Gene", "Func_refGene", "ExonicFunc_refGene", "Gene_ensGene", "avsnp147", "AAChange_ensGene", "AAChange_refGene", "Clinvar", "InterVar_InterVarandEvidence", "Freq_gnomAD_genome_ALL", "Freq_esp6500siv2_all", "Freq_1000g2015aug_all", "CADD_raw", "CADD_phred", "SIFT_score", "GERP++_RS", "phyloP46way_placental" ,"dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "Interpro_domain", "AAChange_knownGene", "rmsk" ,"MetaSVM_score" ,"Freq_gnomAD_genome_POPs", "OMIM" ,"Phenotype_MIM", "OrphaNumber", "Orpha", "Otherinfo"])
    except OSError:
        print ("Could not open/read the input file: " + interVarFinalFileName)

    famid = args.sampleid[3:5]


    hpo_gene_dict = createGeneSyndromeDict(hpo_genes_df)
    hpo_syndrome_dict = createGeneSyndromeDict(hpo_syndromes_df)


    try:
        omim_gene = pd.read_csv(omim_gene, sep='\t', usecols=['Name', 'Score'])
        omim_gene.columns = ['gene', 'syndrome_name']
        omim_gene = omim_gene.loc[omim_gene['syndrome_name']!='.']
        omim_gene = omim_gene.loc[~omim_gene['gene'].isin(hpo_gene_dict.keys())] # This is a list of genes with OMIM disease pheno but not annotated by HPO
        omim_gene = pd.DataFrame({'gene':omim_gene['gene'], 'score':np.zeros(len(omim_gene['gene']))})
    except OSError:
        print("Could not open/read the input file: " + omim_gene)
        sys.exit()


    weights_gene = args.workdir + "/HPO_weight_gene.txt"
    weights_syndrome = args.workdir + "/HPO_weight_syndrome.txt"

    if not os.path.exists(weights_gene):
        cmd = 'grep -v ^# ' + args.database + '/phenotype_to_genes.txt |cut -f-1 | uniq -c |awk \'{print $2, 1/$1}\' > ' + args.workdir + '/HPO_weight_gene.txt'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise Exception(stderr)
    if not os.path.exists(weights_syndrome):
        cmd = 'cat ' + args.database + '/phenotype_annotation.tab | awk -F"\t" \'{print $5}\'| sort | uniq -c | awk \'{print $2, 1/$1}\' > ' + args.workdir + '/HPO_weight_syndrome.txt'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise Exception(stderr)

    weightGeneDict = createWeightDict(weights_gene)
    weightSyndromeDict = createWeightDict(weights_syndrome)

    # Retrieve clinical phenome from the patient
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Analyzing the clinical phenome for ' + args.sampleid + '...')
    clinical_phenome = getClinicalPhenome(args)

    # Get gene sum score
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Generating a clinically relevant primary gene list for ' + args.sampleid + '...')
    gene_score_result_r = calculateGeneSumScore(args, hpo_gene_dict, weightGeneDict, clinical_phenome, omim_gene)

    # Overlap important genes (gene_score_result_r) with all the SNPs and indels
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting SNPs and indels on ' + args.sampleid + '...')
    all_small_variants = smallVariantGeneOverlapCheckInheritance(args, smallVariantFile, interVarFinalFile, gene_score_result_r, famid)


    # If bionano is flagged, check SV from bionano SV calls
    if args.bionano:

        if not ((args.enzyme == 'DLE') or (args.enzyme == 'BspQI')):
            raise Exception("Enzyme flag not set correctly! Either DLE or BspQI must be specified when using -b.")

        # Call all SVs using bionano
        cyto_BN_del, cyto_BN_dup, exon_calls_BN_del, exon_calls_BN_dup = bionanoSV(args, famid, gene_score_result_r, all_small_variants)

    #Make 10x SV calls
    if args.linkedreadSV:
        cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, exon_calls_10x_del, exon_calls_10x_largeSV_del, exon_calls_10x_largeSV_dup = linkedreadSV(args, famid, gene_score_result_r, all_small_variants)

        # Remove control files
        cmd = 'rm ./results/' + args.sampleid + '/10x_del_control.vcf.gz ./results/' + args.sampleid + '/10x_largeSV_control.vcf.gz'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise Exception(stderr)

    # Get differential diagnosis
    if args.bionano or args.linkedreadSV:
        if args.bionano and not args.linkedreadSV:
            cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV = None, None, None
        elif args.linkedreadSV and not args.bionano:
            cyto_BN_del, cyto_BN_dup = None, None

        syndrome_score_result_r = differentialDiangosis(hpo_syndrome_dict, weightSyndromeDict, clinical_phenome, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup, hpo_syndromes_mim_df)

        # Remove control files
        cmd = 'rm ./results/' + args.sampleid + '/bionano_control.smap.gz'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise Exception(stderr)

    if args.bionano and args.linkedreadSV:
        print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Searching for confident large dels and dups on ' + args.sampleid + '...')
        findConfDelDup(args, exon_calls_10x_del, exon_calls_10x_largeSV_del, exon_calls_10x_largeSV_dup, exon_calls_BN_del, exon_calls_BN_dup)

    # Move all the intermediate files to the misc folder
    if not args.singleton:
        cmd = 'mv ./results/' + args.sampleid + '/' + args.sampleid + '_gene_list.txt ./results/' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt ./results/' + args.sampleid + '/' + args.sampleid + '_hpo_inexact.txt ./results/' + args.sampleid + '/' + args.sampleid + '_maternal_inherited_smallVariants.vcf ./results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt ./results/' + args.sampleid + '/' + args.sampleid + '_paternal_inherited_smallVariants.vcf ./results/' + args.sampleid + '/' + args.sampleid + '_smallVariant_candidates.txt ./results/' + args.sampleid + '/' + args.sampleid + '_smallVariants_ALL_candidates.txt ./results/' + args.sampleid + '/' + args.sampleid + '_syndrome_score_result_r.txt ./results/' + args.sampleid + '/' + args.sampleid + '_target.bed ./results/' + args.sampleid + '/' + args.sampleid + '.txt ./results/' + args.sampleid + '/misc/'
    else:
        cmd = 'mv ./results/' + args.sampleid + '/' + args.sampleid + '_gene_list.txt ./results/' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt ./results/' + args.sampleid + '/' + args.sampleid + '_hpo_inexact.txt ./results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt ./results/' + args.sampleid + '/' + args.sampleid + '_smallVariant_candidates.txt ./results/' + args.sampleid + '/' + args.sampleid + '_smallVariants_ALL_candidates.txt ./results/' + args.sampleid + '/' + args.sampleid + '_syndrome_score_result_r.txt ./results/' + args.sampleid + '/' + args.sampleid + '_target.bed ./results/' + args.sampleid + '/' + args.sampleid + '.txt ./results/' + args.sampleid + '/misc/'

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        raise Exception(stderr)

    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Clinical interpretor finished successfully')


    # del_df_name = args.workdir + '/results/' + args.sampleid + "/" + args.sampleid + '_input1.txt'
    # del_df = pd.read_csv(del_df_name, sep='\t')
    # cyto_BN_del_name = args.workdir + '/results/' + args.sampleid + "/" + args.sampleid + '_input2.txt'
    # cyto_BN_del = pd.read_csv(cyto_BN_del_name, sep='\t')
    # #cyto_BN_del = cyto_BN_del.drop(columns=['CHROM_b', 'POS_b', 'ID_b', 'REF_b', 'ALT_1_b','ALT_2_b','ALT_3_b', 'QUAL_b','FILTER_PASS_b', 'END_b','SVLEN_b','Start_b', 'End_b'])
    # cols = []
    # delDupSyndromeSVOverlap(del_df, cyto_BN_del, cols)


if __name__=="__main__":
    main()


