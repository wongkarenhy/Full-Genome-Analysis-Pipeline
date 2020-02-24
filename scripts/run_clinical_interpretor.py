#!/usr/bin/env python3.6

import pandas as pd
from collections import defaultdict, Counter
import argparse
import sys
import os
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
    for var, hpo in database_df.itertuples(index=False): # var can either be gene or syndrome
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



def calculateGeneSumScore(args, hpo_gene_dict, weightDict, clinical_phenome):

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
            hpo_pheno = set(hpo_gene_dict[query]) # To get the phenotypic features for a given gene
            overlap = hpo_pheno.intersection(clinical_phenome) # overlap all the phenotypic features with the clinical phenomes

            for term in overlap:
                gene_sum_score += weightDict[term]

            gene_score_result = gene_score_result.append({'gene':query, 'score':gene_sum_score}, ignore_index=True)

        gene_score_result_r = gene_score_result.iloc[::-1]
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
    #df = df[['Chr', 'Start', 'End', 'Clinvar', 'InterVar_InterVarandEvidence']]
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
    df = df.sort_values(by='Patho_score', ascending=False)
    df = df.drop(['Clinvar_idx', 'InterVar_idx', 'Clinvar_score', 'InterVar_score', 'Patho_score'], axis=1)

    return df



def smallVariantGeneOverlapCheckInheritance(args, smallVariantFile, interVarFinalFile, gene_score_result_r, famid):

    # Overlap gene_score_result_r with small variants genes found in the proband
    gene_score_result_r = gene_score_result_r[gene_score_result_r.gene.isin(smallVariantFile.gene)]

    #gene_score_result_r = normalizeRawScore(args, gene_score_result_r, 'gene')

    # Subset the intervar files further to store entries relevant to these set of genes
    filtered_intervar = pd.merge(interVarFinalFile, gene_score_result_r, left_on='Ref_Gene', right_on='gene',how='inner')
    filtered_intervar = filtered_intervar.sort_values(by='score', ascending=False)


    # Create a bed file and write it out
    pd.DataFrame(filtered_intervar).to_csv(
        './results/' + args.sampleid + "/" + args.sampleid + '_smallVariant_candidates.txt', index=False, sep='\t',
        header=False)  # Write out a subset of the variant first
    filtered_intervar_bed = filtered_intervar[['Chr', 'Start', 'End']]
    filtered_intervar_bed.loc[:,'Chr'] = 'chr' + filtered_intervar_bed.loc[:,'Chr'].astype(str)
    filtered_intervar_bed.loc[:,'Start'] -= 1
    pd.DataFrame(filtered_intervar_bed).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_target.bed', index=False, sep='\t', header=False)

    # Get overlapping variants from the parents so we know which variants are inherited
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Comparing small variants (SNPs/indels) inheritance')
    cmd = "bcftools view -R ./results/" + args.sampleid + "/" + args.sampleid + "_target.bed /media/KwokRaid04/CIAPM/CIAPM_longranger/BC0" + famid + "01_longranger/outs/phased_variants.vcf.gz > ./results/" + args.sampleid + "/" + args.sampleid + "_paternal_inherited_smallVariants.vcf"
    os.system(cmd)
    cmd = "bcftools view -R ./results/" + args.sampleid + "/" + args.sampleid + "_target.bed /media/KwokRaid04/CIAPM/CIAPM_longranger/BC0" + famid + "02_longranger/outs/phased_variants.vcf.gz > ./results/" + args.sampleid + "/" + args.sampleid + "_maternal_inherited_smallVariants.vcf"
    os.system(cmd)

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
    #dominant = filtered_intervar[((filtered_intervar['paternal'] == 1) & (filtered_intervar['maternal'] == 0)) | ((filtered_intervar['maternal'] == 1) & (filtered_intervar['paternal'] == 0))]
    ## De novo
    denovo = filtered_intervar[(filtered_intervar['paternal'] == 0) & (filtered_intervar['maternal'] == 0)]
    #Compound het
    filtered_intervar_compoundhet = filtered_intervar[(filtered_intervar['Otherinfo'] == 'het')]
    filtered_intervar_compoundhet = filtered_intervar_compoundhet[(filtered_intervar_compoundhet['maternal'] != 2) & (filtered_intervar_compoundhet['paternal'] != 2) & ((filtered_intervar_compoundhet['paternal'] == 1) & (filtered_intervar_compoundhet['maternal'] == 0)) | ((filtered_intervar_compoundhet['maternal'] == 1) & (filtered_intervar_compoundhet['paternal'] == 0))]
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
    pd.DataFrame(recessive).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_smallVariants_recessive_candidates.txt', index=False, sep='\t', header=True)
    # Dominant
    #pd.DataFrame(dominant).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_dominant_smallVariants_candidates.txt', index=False, sep='\t', header=True)
    # De novo
    pd.DataFrame(denovo).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_smallVariants_denovo_candidates.txt', index=False, sep='\t', header=True)
    # Compound het
    pd.DataFrame(compoundhet).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_smallVariants_compoundhet_candidates.txt', index=False, sep='\t', header=True)
    # All
    pd.DataFrame(filtered_intervar).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_smallVariants_ALL_candidates.txt', index=False, sep='\t', header=True)

    return filtered_intervar



def differentialDiangosis(hpo_syndrome_dict, weightSyndromeDict, clinical_phenome, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup):

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
    delDupSyndrome(syndrome_score_result_r, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup)
    #delDupSyndrome(syndrome_score_result_r, args)

    return(syndrome_score_result_r)




def findGenomicLocation(cytoband_key, cytobandDict):

    keys = [key for key in cytobandDict if key.startswith(cytoband_key)]
    genomic_coords_list = []
    for key in keys:
        genomic_coords_list.append(cytobandDict[key].split('-'))

    genomic_coords_list = list(chain.from_iterable(genomic_coords_list))
    min_coords = min(genomic_coords_list)
    max_coords = max(genomic_coords_list)

    genomic_range = str(min_coords) + '-' + str(max_coords)

    return genomic_range





def parseSyndromeNameToCytoband(df, cytobandDict):

    regex = r'([0-9XY]{1,2}[PQ]{1}[\w\\.\\-]{1,10}[\s$])'

    for index, row in df.iterrows():
        m = re.search(regex, str(row))
        if m is not None:
            df.loc[index, 'cytoband'] = m.group(1)

    df.dropna(subset=['cytoband'], inplace=True)
    df['cytoband'] = df['cytoband'].str.lower()
    df['cytoband'] = df['cytoband'].str.replace('x', 'X')
    df[['Chromosome', 'discard']] = df.cytoband.str.split('p|q', 1, expand=True)
    df = df.drop('discard', axis=1)
    df[['cytoband_start', 'cytoband_stop']] = df.cytoband.str.split('-', expand=True)
    df['arm'] = np.where(df['cytoband_start'].str.contains('p'), 'p', 'q')

    for idx, row in df.iterrows():

        cytoband_start_key = row['cytoband_start'].replace(" ","")
        if cytoband_start_key in cytobandDict:
            coords_start = cytobandDict[cytoband_start_key]

        else:
            genomic_range = findGenomicLocation(cytoband_start_key, cytobandDict)
            coords_start = genomic_range


        if row['cytoband_stop'] is not None: # Fix cytoband_stop column for quick cytobandDict lookup
            print(type(row['cytoband_stop']))
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




def delDupSyndrome(syndrome_score_result_r, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup):

    # Create cytoband <-> genomic coordinates dict
    cytobandDict = createCytobandDict(args)

    del_cond = syndrome_score_result_r['syndrome'].str.contains('DELETION')
    dup_cond = syndrome_score_result_r['syndrome'].str.contains('DUPLICATION')

    del_df = syndrome_score_result_r[del_cond]
    dup_df = syndrome_score_result_r[dup_cond]

    del_df = parseSyndromeNameToCytoband(del_df, cytobandDict)
    dup_df = parseSyndromeNameToCytoband(dup_df, cytobandDict)

    if args.bionano:

        # Overlap with del/dup syndromes
        if cyto_BN_dup is not None: # It can be None because old Bionano pipeline doesn't call duplications...
            overlap_dup_BN = delDupSyndromeSVOverlap(dup_df, cyto_BN_dup)
            overlap_dup_BN.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_duplication_syndrome.txt', sep='\t', index=False)
        else:
            pd.DataFrame().to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_duplication_syndrome.txt', sep='\t', index=False)

        overlap_del_BN = delDupSyndromeSVOverlap(del_df, cyto_BN_del)
        overlap_del_BN.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_deletion_syndrome.txt', sep='\t', index=False)

    if args.linkedreadSV:
        overlap_dup_largeSV_10x = delDupSyndromeSVOverlap(dup_df, cyto_10x_dup_largeSV)
        overlap_dup_largeSV_10x.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_duplication_largeSV_syndrome.txt', sep='\t', index=False)

        overlap_del_largeSV_10x = delDupSyndromeSVOverlap(del_df, cyto_10x_del_largeSV)
        overlap_del_largeSV_10x.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_deletion_largeSV_syndrome.txt', sep='\t', index=False)

        overlap_del_10x = delDupSyndromeSVOverlap(del_df, cyto_10x_del)
        overlap_del_10x.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_deletion_syndrome.txt', sep='\t', index=False)

    #if args.linkedreadSV and args.bionano:
        # syndrome appearing in both 10x and bionano --> confident set





def delDupSyndromeSVOverlap(del_df, cyto_BN_del):
    overlap_del_BN = PyRanges(del_df).overlap(PyRanges(cyto_BN_del))
    if not overlap_del_BN.df.empty:
        overlap_del_BN = overlap_del_BN.df.sort_values(by='score', ascending=False)
        # Rearrange the column
        cols = ['Chromosome', 'Start', 'End', 'cytoband', 'cytoband_start', 'cytoband_stop', 'arm', 'syndrome', 'score', 'normalized_score']
        overlap_del_BN = overlap_del_BN[cols]
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

    return full_paths





def bionanoSV(args, famid, gene_score_result_r, all_small_variants):

    # Generate controls files (1KGP BN samples + CIAPM parents (excluding  parents of the proband of interest)
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Generating bionano control file...')
    control_files_path = [args.workdir + "/bionano_sv/controls/DLE", args.workdir + "/bionano_sv/controls/BspQI", args.workdir + "/bionano_sv/cases/DLE", args.workdir + "/bionano_sv/cases/BspQI"]
    full_paths = compileControlFiles(control_files_path, famid)

    ## Write an empty file
    with open(args.workdir + "/results/" + args.sampleid + "/bionano_control.smap.gz", 'w'):  # So it will overwrite the old file
        pass

    for path in full_paths:
        cmd = "cat " + path + "/exp_refineFinal1_merged_filter.smap | gzip >> " + args.workdir + "/results/" + args.sampleid + "/bionano_control.smap.gz"
        os.system(cmd)

    # Create a BN arg object
    BN_args = Namespace(sampleID = args.sampleid,
                        samplepath = args.workdir + "/bionano_sv/cases/" + args.enzyme + "/" + args.sampleid + "/exp_refineFinal1_merged_filter.smap",
                        fpath = args.workdir + "/bionano_sv/cases/" + args.enzyme + "/BC0" + famid + "01/exp_refineFinal1_merged_filter.smap",
                        mpath = args.workdir + "/bionano_sv/cases/" + args.enzyme + "/BC0" + famid + "02/exp_refineFinal1_merged_filter.smap",
                        referencepath = args.workdir + "/results/" + args.sampleid + "/bionano_control.smap.gz",
                        outputdirectory = args.workdir + '/results/' + args.sampleid,
                        exons = args.workdir + '/annotatedexonsphenotypes.bed',
                        genelist = gene_score_result_r)


    # Call bionano translocation
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano translocations on ' + args.sampleid + '...')
    BN_translocation(BN_args)

    # Call bionano deletion
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano deletions on ' + args.sampleid + '...')
    cyto_BN_del, exon_calls_BN_del = BN_deletion(BN_args)

    # Call bionano insertion
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano insertions on ' + args.sampleid + '...')
    BN_insertion(BN_args)

    # Call bionano duplications
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano duplications on ' + args.sampleid + '...')
    cyto_BN_dup, exon_calls_BN_dup = BN_duplication(BN_args)

    # Call bionano inversions
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting bionano inversions on ' + args.sampleid + '...')
    BN_inversion(BN_args)

    # Check potential compoundhets with SNPs and indels
    #all_small_variants_exons = all_small_variants['Ref_Gene']
    BN_exons = pd.concat([exon_calls_BN_del, exon_calls_BN_dup])
    if BN_exons.empty:
        pd.DataFrame().to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_SV_SNPsIndels_compounthet_candidates.txt', sep='\t', index=False)
    else:
        #BN_exons = BN_exons[BN_exons['gene'].isin(all_small_variants_exons)]
        BN_exons = pd.merge(BN_exons, all_small_variants, left_on='gene', right_on='Ref_Gene', how='inner')
        BN_exons.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_Bionano_SV_SNPsIndels_compounthet_candidates.txt', sep='\t', index=False)

    return cyto_BN_del, cyto_BN_dup




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
        os.system(cmd)


    # Need to generate another reference file for large SVs
    with open(args.workdir + "/results/" + args.sampleid + "/10x_largeSV_control.vcf.gz",'w'):  # So it will overwrite the old file
        pass

    for path in full_paths:
        cmd = "zcat " + path + "/large_svs.vcf.gz | gzip >> " + args.workdir + "/results/" + args.sampleid + "/10x_largeSV_control.vcf.gz"
        os.system(cmd)

    tenx_args_del = Namespace(sampleID = args.sampleid,
                        samplepath = args.workdir + "/linkedRead_sv/cases/" + args.sampleid + "/dels.vcf.gz",
                        fpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "01/dels.vcf.gz",
                        mpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "02/dels.vcf.gz",
                        referencepath = args.workdir + "/results/" + args.sampleid + "/10x_del_control.vcf.gz",
                        outputdirectory = args.workdir + '/results/' + args.sampleid,
                        exons = args.workdir + '/annotatedexonsphenotypes.bed',
                        genelist = gene_score_result_r)

    tenx_args_largeSV = Namespace(sampleID = args.sampleid,
                        samplepath = args.workdir + "/linkedRead_sv/cases/" + args.sampleid + "/large_svs.vcf.gz",
                        fpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "01/large_svs.vcf.gz",
                        mpath = args.workdir + "/linkedRead_sv/cases/BC0" + famid + "02/large_svs.vcf.gz",
                        referencepath = args.workdir + "/results/" + args.sampleid + "/10x_largeSV_control.vcf.gz",
                        outputdirectory = args.workdir + '/results/' + args.sampleid,
                        exons = args.workdir + '/annotatedexonsphenotypes.bed',
                        genelist = gene_score_result_r)

    # Call medium size deletions
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads medium deletions on ' + args.sampleid + '...')
    cyto_10x_del, exon_calls_10x_del = tenxdeletions(tenx_args_del)

    # Call large deletions
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large deletions on ' + args.sampleid + '...')
    cyto_10x_del_largeSV, exon_calls_10x_largeSV_del = tenxlargesvdeletions(tenx_args_largeSV)

    # Call large duplications
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large duplications on ' + args.sampleid + '...')
    cyto_10x_dup_largeSV, exon_calls_10x_largeSV_dup = tenxlargesvduplications(tenx_args_largeSV)

    # Call large inversions
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting linked-reads large inversions on ' + args.sampleid + '...')
    tenxlargesvinversions(tenx_args_largeSV)

    # Check potential compoundhets with SNPs and indels
    tenx_exons = pd.concat([exon_calls_10x_del, exon_calls_10x_largeSV_del, exon_calls_10x_largeSV_dup])
    if tenx_exons.empty:
        pd.DataFrame().to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_SV_SNPsIndels_compounthet_candidates.txt', sep='\t', index=False)
    else:
        #BN_exons = BN_exons[BN_exons['gene'].isin(all_small_variants_exons)]
        tenx_exons = pd.merge(tenx_exons, all_small_variants, left_on='gene', right_on='Ref_Gene', how='inner')
        tenx_exons.to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_10x_SV_SNPsIndels_compounthet_candidates.txt', sep='\t', index=False)

    return cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV



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
    args = parser.parse_args()

    # Change work dir
    os.chdir(args.workdir)

    # Define variables
    ## Read the database files
    hpo_genes = args.database + "/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt"
    hpo_syndromes = args.database + "/phenotype_annotation_hpoteam.tab"
    smallVariantFileName = args.intervar + "/example/" + args.sampleid + "_smallVariant_geneList.txt"
    interVarFinalFileName = args.intervar + "/example/" + args.sampleid +  ".hg38_multianno.txt.intervar.FINAL"

    try:
        hpo_genes_df = pd.read_csv(hpo_genes, sep='\t', usecols=[1, 3], names=["gene_name", "HPO_id"], comment='#')
    except OSError:
        print("Could not open/read the input file: " + hpo_genes)
        sys.exit()

    try:
        hpo_syndromes_df = pd.read_csv(hpo_syndromes, sep='\t', usecols=[2, 4], names=["syndrome_name", "HPO_id"])
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


    weights_gene = "./HPO_weight_gene.txt"
    weights_syndrome = "./HPO_weight_syndrome.txt"
    famid = args.sampleid[3:5]

    hpo_gene_dict = createGeneSyndromeDict(hpo_genes_df)
    hpo_syndrome_dict = createGeneSyndromeDict(hpo_syndromes_df)

    weightGeneDict = createWeightDict(weights_gene)
    weightSyndromeDict = createWeightDict(weights_syndrome)

    # Retrieve clinical phenome from the patient
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Analyzing the clinical phenome for ' + args.sampleid + '...')
    clinical_phenome = getClinicalPhenome(args)

    # Get gene sume score
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Generating a clinically relevant primary gene list for ' + args.sampleid + '...')
    gene_score_result_r = calculateGeneSumScore(args, hpo_gene_dict, weightGeneDict, clinical_phenome)

    # Overlap important genes (gene_score_result_r) with all the SNPs and indels
    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Detecting SNPs and indels on ' + args.sampleid + '...')
    all_small_variants = smallVariantGeneOverlapCheckInheritance(args, smallVariantFile, interVarFinalFile, gene_score_result_r, famid)


    # If bionano is flagged, check SV from bionano SV calls
    if args.bionano:

        if not ((args.enzyme == 'DLE') or (args.enzyme == 'BspQI')):
            raise Exception("Enzyme flag not set correctly! Either DLE or BspQI must be specified when using -b.")

        # Call all SVs using bionano
        cyto_BN_del, cyto_BN_dup = bionanoSV(args, famid, gene_score_result_r, all_small_variants)

    #Make 10x SV calls
    if args.linkedreadSV:
        cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV = linkedreadSV(args, famid, gene_score_result_r, all_small_variants)

        # Remove control files
        cmd = 'rm ./results/' + args.sampleid + '/10x_del_control.vcf.gz ./results/' + args.sampleid + '/10x_largeSV_control.vcf.gz'
        os.system(cmd)

    # Get differential diagnosis
    if args.bionano or args.linkedreadSV:
        if args.bionano and not args.linkedreadSV:
            cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV = None, None, None
        elif args.linkedreadSV and not args.bionano:
            cyto_BN_del, cyto_BN_dup = None, None

        syndrome_score_result_r = differentialDiangosis(hpo_syndrome_dict, weightSyndromeDict, clinical_phenome, args, cyto_10x_del, cyto_10x_del_largeSV, cyto_10x_dup_largeSV, cyto_BN_del, cyto_BN_dup)

        # Remove control files
        cmd = 'rm ./results/' + args.sampleid + '/bionano_control.smap.gz'
        os.system(cmd)

    print('[run_clinical_interpretor.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Pipeline finished successfully')


if __name__=="__main__":
    main()






