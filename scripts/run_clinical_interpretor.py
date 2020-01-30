#!/user/bin/python3
import pandas as pd
from collections import defaultdict
import argparse
import sys

def createGeneSyndromeDict(database_df):

    dict = defaultdict(list)
    for var, hpo in database_df.itertuples(index=False): # var can either be gene or syndrome
        dict[var].append(hpo)

    return(dict)



def createWeightDict(weights):

    try:
        w_df = pd.read_csv(weights, sep = ' ', names=["HPO_id", "weight"], comment = '#')
    except OSError:
        print("Count not open/read the input file:")
        sys.exit()

    weightDict = dict(zip(w_df.HPO_id, w_df.weight))

    return(weightDict)



def getClinicalPhenome(args):
    # Get the clinical phenome and store as a set
    try:
        clinical_phenome = set(open("./genelist/" + args.sampleid + "_hpo_inexact.txt").read().splitlines())
    except OSError:
        print("Count not open/read the input file:")
        sys.exit()

    return(clinical_phenome)



def calculateGeneSumScore(args, hpo_gene_dict, weightDict, clinical_phenome):

    # Go through genes in genelist found in the patients
    try:
        genes = open("./genelist/" + args.sampleid + "_gene_list.txt", 'r')
    except OSError:
        print("Count not open/read the input file:")
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

    return(syndrome_score_result_r)



def normalizeRawScore(args, raw_score, mode):

    # Normalize all the scores to 1-100
    max_score = max(raw_score['score'])
    raw_score['normalized_score'] = raw_score['score']/max_score * 100

    pd.DataFrame(raw_score).to_csv('./genelist/' + args.sampleid + '_' + mode + '_score_result.txt', index=False, sep='\t', header=True)



def main():

    # Parse argument
    parser = argparse.ArgumentParser(description="This software ranks genes based on the clinical phenome.")
    parser.add_argument("-s", "--sampleid",
                        help="Sample ID",
                        dest="sampleid", type=str)
    args = parser.parse_args()

    # Define variables
    ## Read the database files
    hpo_genes = "../../../../../media/KwokRaid02/karen/database/human_pheno_ontology/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt"
    try:
        hpo_genes_df = pd.read_csv(hpo_genes, sep='\t', usecols=[1, 3], names=["gene_name", "HPO_id"], comment='#')
    except OSError:
        print("Count not open/read the input file:")
        sys.exit()

    hpo_syndromes = "../../../../../media/KwokRaid02/karen/database/human_pheno_ontology/phenotype_annotation_hpoteam.tab"
    try:
        hpo_syndromes_df = pd.read_csv(hpo_syndromes, sep='\t', usecols=[2, 4], names=["syndrome_name", "HPO_id"])
        #hpo_syndromes_df = pd.read_csv(hpo_syndromes, sep='\t', usecols=[1,2,3,4], names=["col1", "col2", "col3", "col4"])
    except OSError:
        print("Count not open/read the input file:")
        sys.exit()

    #print(hpo_syndromes_df)

    weights_gene = "./HPO_weight_gene.txt"
    weights_syndrome = "./HPO_weight_syndrome.txt"

    hpo_gene_dict = createGeneSyndromeDict(hpo_genes_df)
    hpo_syndrome_dict = createGeneSyndromeDict(hpo_syndromes_df)

    weightGeneDict = createWeightDict(weights_gene)
    weightSyndromeDict = createWeightDict(weights_syndrome)

    # Retrieve clinical phenome from the patient
    clinical_phenome = getClinicalPhenome(args)

    # Get gene sume score and differential diagnosis
    gene_score_result_r = calculateGeneSumScore(args, hpo_gene_dict, weightGeneDict, clinical_phenome)

    # Get differential diagnosis
    syndrome_score_result_r = differentialDiangosis(hpo_syndrome_dict, weightSyndromeDict, clinical_phenome)

    normalizeRawScore(args, gene_score_result_r, 'gene')
    normalizeRawScore(args, syndrome_score_result_r, 'syndrome')


if __name__=="__main__":
    main()








