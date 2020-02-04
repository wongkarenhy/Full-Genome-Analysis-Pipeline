#!/user/bin/python3
import obonet
import networkx
import argparse
import pandas as pd
import os

def queryHPO(args):

    relatives = set()

    # Read the hpo file in obo format
    try:
        hpo = obonet.read_obo(args.database + "/hp.obo")
    except OSError:
        print("Count not open/read the obo file:" + args.database + "/hp.obo")
        sys.exit()

    # Read the file containing query HPO terms
    try:
        hpo_exact = open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt', 'r')
    except OSError:
        print("Count not open/read the input file:", args.workdir + '/results' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt')
        sys.exit()

    # Read the file containing manual HPO terms
    try:
        hpo_manual = open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt', 'r')
    except OSError:
        print("Count not open/read the input file:", args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt')
        sys.exit()

    with hpo_manual:

        hpo_manual_file = hpo_manual.read().splitlines()

        for manual_term in hpo_manual_file:

            relatives.add(manual_term)

    with hpo_exact:

        hpo_exact_file = hpo_exact.read().splitlines()

        for query in hpo_exact_file:

            # Get the parents and children terms
            # Relatives include self
            relatives.add(query)
            relatives = relatives.union(networkx.ancestors(hpo, query))
            relatives = relatives.union(networkx.descendants(hpo, query))

        with open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_inexact.txt', 'a') as out:

            for terms in relatives:

                out.write(terms + '\n')

    return(relatives)

def getGeneList(args, relatives):

    geneList = set()

    try:
        # Read the gene phenotype file into a pandas dataframe
        df = pd.read_csv(args.database + "/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt", sep = '\t', names=["HPO_id", "HPO_name", "gene_id", "gene_name"], comment = '#')
    except OSError:
        print("Count not open/read the gene-phenotype file:", args.database + "/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt")
        sys.exit()

    for hpo in relatives:

        geneList = geneList.union(set(df[df.HPO_id == hpo]['gene_name'].tolist()))

    with open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_gene_list.txt', 'a') as out:

        for gene in geneList:

            out.write(gene + '\n')


def main():

    parser=argparse.ArgumentParser(description="This software extracts inexact HPO terms and generate a gene list.")
    parser.add_argument("-s","--sampleid",help="Sample ID." ,dest="sampleid", type=str, required = True)
    parser.add_argument("-w","--workdir",help="This is the base work directory (folder/directory).",dest="workdir",type=str, required = True)
    parser.add_argument("-d", "--database", help="Path to the HPO database", dest="database", type=str, required = True)
    # parser.add_argument("-o", "--ontology", \
    #                     help="Ontology file in obo format. Default is /media/KwokRaid02/karen/database/human_pheno_ontology/hp.obo", \
    #                     dest="ontology", type=str, default = "/media/KwokRaid02/karen/database/human_pheno_ontology/hp.obo")
    # parser.add_argument("-g", "--gene", \
    #                     help="Gene phenotype relationship table. Default is /media/KwokRaid02/karen/database/human_pheno_ontology/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt", \
    #                     dest="gene", type=str, default = "/media/KwokRaid02/karen/database/human_pheno_ontology/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt")
    args=parser.parse_args()


    # run query
    relatives = queryHPO(args)
    getGeneList(args, relatives)

if __name__=="__main__":
    main()




