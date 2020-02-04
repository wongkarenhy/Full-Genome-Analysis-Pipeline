#!/user/bin/python3
import pandas as pd
from collections import defaultdict, Counter
import argparse
import sys
import os
import re
import numpy as np

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

    return(gene_score_result_r)



def getParentsGeno(filtered_intervar, inheritance_mode, ov_allele):

    # Create two new columns and initialize to 0
    filtered_intervar[inheritance_mode] = 0
    filtered_intervar = filtered_intervar.reset_index(drop=True)

    for idx, row in enumerate(filtered_intervar.itertuples(index=False)):
        if int(getattr(row, 'Start')) in set(ov_allele['Start']):
            parents_geno = ov_allele.loc[ov_allele['Start']==getattr(row,'Start'),'geno'].item()
            filtered_intervar.set_value(idx, inheritance_mode, parents_geno)

    return(filtered_intervar)



def smallVariantGeneOverlapCheckInheritance(args, smallVariantFile, interVarFinalFile, gene_score_result_r):

    # Overlap gene_score_result_r with small variants genes found in the proband
    gene_score_result_r = gene_score_result_r[gene_score_result_r.gene.isin(smallVariantFile.gene)]

    gene_score_result_r = normalizeRawScore(args, gene_score_result_r, 'gene')

    # Subset the intervar files further to store entries relevant to these set of genes
    ### BC00103.hg38_multianno.txt.intervar.FINAL
    # filtered_intervar = interVarFinalFile[interVarFinalFile.Ref_Gene.isin(gene_score_result_r.gene)]
    filtered_intervar = pd.merge(interVarFinalFile, gene_score_result_r, left_on='Ref_Gene', right_on='gene',how='inner')
    filtered_intervar = filtered_intervar.sort_values(by='score', ascending=False)

    # Create a bed file and write it out
    pd.DataFrame(filtered_intervar).to_csv(
        './results/' + args.sampleid + "/" + args.sampleid + '_smallVariantCandidates.txt', index=False, sep='\t',
        header=False)  # Write out a subset of the variant first
    filtered_intervar_bed = filtered_intervar[['Chr', 'Start', 'End']]
    filtered_intervar_bed.loc[:,'Chr'] = 'chr' + filtered_intervar_bed.loc[:,'Chr'].astype(str)
    filtered_intervar_bed.loc[:,'Start'] -= 1
    pd.DataFrame(filtered_intervar_bed).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_target.bed', index=False, sep='\t', header=False)

    # Get overlapping variants from the parents so we know which variants are inherited
    famid = args.sampleid[3:5]
    print("Comparing small variants (SNPs/indels) inheritance")
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
        col_list = ['paternal', 'maternal']
        res = df[col_list].sum(axis=0)
        if (res[0] == 0) or (res[1] == 0):
            discard.append(gene)

    compoundhet = compoundhet[~compoundhet['Ref_Gene'].isin(discard)]


    # Print all the variants according to inheritance mode
    # Recessive
    pd.DataFrame(recessive).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_recessive_smallVariants_candidates.txt', index=False, sep='\t', header=True)
    # Dominant
    #pd.DataFrame(dominant).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_dominant_smallVariants_candidates.txt', index=False, sep='\t', header=True)
    # De novo
    pd.DataFrame(denovo).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_denovo_smallVariants_candidates.txt', index=False, sep='\t', header=True)
    # Compound het
    pd.DataFrame(compoundhet).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_compoundhet_smallVariants_candidates.txt', index=False, sep='\t', header=True)
    # All
    pd.DataFrame(filtered_intervar).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_ALL_smallVariants_candidates.txt', index=False, sep='\t', header=True)



def differentialDiangosis(hpo_syndrome_dict, weightSyndromeDict, clinical_phenome):

    syndrome_score_result = pd.DataFrame(columns=['syndrome', 'score'])


    # Check every syndrome and its overlapping hpo terms
    for syndrome in hpo_syndrome_dict:
        #print(syndrome)
        hpo_terms = set(hpo_syndrome_dict[syndrome])

        score = 0
        for term in hpo_terms:

            if term in clinical_phenome:
                score += weightSyndromeDict[term]

        if score != 0:
            syndrome_score_result = syndrome_score_result.append({'syndrome': syndrome, 'score': score}, ignore_index=True)

    syndrome_score_result_r = syndrome_score_result.sort_values(by='score', ascending=False)
    syndrome_score_result_r['syndrome'] = syndrome_score_result_r['syndrome'].str.upper()
    #print(syndrome_score_result_r)
    # Specifically look for deletion/duplication syndrome
    delDupSyndrome(syndrome_score_result_r)

    return(syndrome_score_result_r)




def parseSyndromeNameToCytoband(df):

    regex = r'([0-9XY]{1,2}[PQ]{1}[\w\\.\\-]{1,10}[\s$])'

    for index, row in df.iterrows():
        m = re.search(regex, str(row))
        if m is not None:
            df.ix[index, 'cytoband'] = m.group(1)

    return(df)



def delDupSyndrome(syndrome_score_result_r):

    del_cond = syndrome_score_result_r['syndrome'].str.contains('DELETION')
    dup_cond = syndrome_score_result_r['syndrome'].str.contains('DUPLICATION')

    del_df = syndrome_score_result_r[del_cond]
    dup_df = syndrome_score_result_r[dup_cond]

    del_df = parseSyndromeNameToCytoband(del_df)
    dup_df = parseSyndromeNameToCytoband(dup_df)

    # Overlap with Michelle's Del/Dup datasets



def normalizeRawScore(args, raw_score, mode):

    # Normalize all the scores to 1-100
    max_score = max(raw_score['score'])
    raw_score.loc[:,'normalized_score'] = raw_score.loc[:,'score']/max_score * 100

    #pd.DataFrame(raw_score).to_csv('./results/' + args.sampleid + "/" + args.sampleid + '_' + mode + '_score_result.txt', index=False, sep='\t', header=True)

    return(raw_score)


def main():

    # Parse argument
    parser = argparse.ArgumentParser(description="This software ranks genes based on the clinical phenome.")
    parser.add_argument("-s", "--sampleid",help="Sample ID",dest="sampleid", type=str, required = True)
    parser.add_argument("-w", "--workdir", help="This is the base work directory.", dest="workdir", type=str, required = True)
    parser.add_argument("-d", "--database", help="Path to HPO database", dest="database", type=str, required = True)
    parser.add_argument("-i", "--intervar", help="Path to InterVar output folder", dest="intervar", type=str, required = True)
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

    hpo_gene_dict = createGeneSyndromeDict(hpo_genes_df)
    hpo_syndrome_dict = createGeneSyndromeDict(hpo_syndromes_df)

    weightGeneDict = createWeightDict(weights_gene)
    weightSyndromeDict = createWeightDict(weights_syndrome)

    # Retrieve clinical phenome from the patient
    clinical_phenome = getClinicalPhenome(args)

    # Get gene sume score
    # Overlap the gene list (gene_score_result_r) with the snv, indel list generated as part of Intervar
    gene_score_result_r = calculateGeneSumScore(args, hpo_gene_dict, weightGeneDict, clinical_phenome)

    # Overlap important genes (gene_score_result_r) with all the SNPs and indels
    smallVariantGeneOverlapCheckInheritance(args, smallVariantFile, interVarFinalFile, gene_score_result_r)

    # Get differential diagnosis
    syndrome_score_result_r = differentialDiangosis(hpo_syndrome_dict, weightSyndromeDict, clinical_phenome)


    #normalizeRawScore(args, syndrome_score_result_r, 'syndrome')



if __name__=="__main__":
    main()








