#!/usr/bin/env python3.6
import obonet
import networkx
import argparse
import pandas as pd
import os
import subprocess
from datetime import datetime


def queryHPO(args):

    relatives = set()

    # Read the hpo file in obo format
    try:
        hpo = obonet.read_obo(args.workdir + "/human_pheno_ontology_b1270/hp.obo")
    except OSError:
        print("Count not open/read the obo file:" + args.workdir + "/human_pheno_ontology_b1270/hp.obo")
        sys.exit()

    # Read the file containing manual HPO terms
    try:
        #hpo_manual = open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt', 'r')
        hpo_manual = pd.read_csv(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt', sep='\t', names = ["HPO_id"])
    except OSError:
        print("Count not open/read the input file:", args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt')
        sys.exit()


    if args.manual:
        all_hpo = set(hpo_manual['HPO_id'])
    else:
        # Read the file containing query HPO terms
        try:
            #hpo_exact = open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt', 'r')
            hpo_exact = pd.read_csv(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt', sep='\t', names = ["HPO_id"])
            all_hpo = set(hpo_exact['HPO_id'].append(hpo_manual['HPO_id']))
        except OSError:
            print("Count not open/read the input file:", args.workdir + '/results' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt')
            sys.exit()


    # with hpo_manual:
    #
    #     hpo_manual_file = hpo_manual.read().splitlines()
    #
    #     for manual_term in hpo_manual_file:
    #
    #         relatives.add(manual_term)

    # with hpo_exact:
    #
    #     hpo_exact_file = hpo_exact.read().splitlines()

    for query in all_hpo:

        # Get the parents and children terms
        # Relatives include self
        relatives.add(query)
        relatives = relatives.union(networkx.ancestors(hpo, query))
        relatives = relatives.union(networkx.descendants(hpo, query))

    with open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_inexact.txt', 'w') as out:

        for terms in relatives:

            out.write(terms + '\n')

    return(relatives)




def getGeneList(args, relatives):

    geneList = set()

    try:
        # Read the gene phenotype file into a pandas dataframe
        df = pd.read_csv(args.workdir + "/human_pheno_ontology_b1270/phenotype_to_genes.txt", sep = '\t', usecols=[*range(0,4)], names=["HPO_id", "HPO_name", "gene_id", "gene_name"], comment = '#')
    except OSError:
        print("Count not open/read the gene-phenotype file:", args.workdir + "/human_pheno_ontology_b1270/phenotype_to_genes.txt")
        sys.exit()

    for hpo in relatives:

        geneList = geneList.union(set(df[df.HPO_id == hpo]['gene_name'].tolist()))

    with open(args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_gene_list.txt', 'w') as out:

        for gene in geneList:

            out.write(gene + '\n')




def run_clinphen(args):

    # Extract text from clinical notes
    #cmd1 = ['cat' , "{}".format(args.json) , '|','jq-linux64', '.features[].notes','|','grep', '-v', 'null', '|', 'sed', '-e', "'s/$/\./g'", '>',  "{}".format(args.workdir),  '/results/',  "{}".format(args.sampleid),  '/',  "{}".format(args.sampleid), '_tmp.txt']

    cmd1 = 'cat ' + args.json + " | jq-linux64 '.features[].notes' | grep -v null | sed -e 's/$/\./g' > " + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_tmp.txt'
    cmd2 = 'cat ' + args.json + " | jq-linux64 '.notes.medical_history' >> " + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_tmp.txt'
    cmd3 = 'cat ' + args.json + " | jq-linux64 '.features[].id' | tr -d '\"' > " + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_manual.txt'

    # clean up the file
    cmd4 = 'cat ' + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + "_tmp.txt | tr -d '[]\"\n' | tr '\' '!' | sed 's/!r!n/\./g' > " + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '.txt'

    # extract HPO terms
    cmd5 = 'clinphen ' + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + ".txt | awk -F'\t' '{print $1}' | tail -n +2 > " + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_hpo_exact.txt'

    # rm tmp files
    cmd6 = 'rm ' + args.workdir + '/results/' + args.sampleid + '/' + args.sampleid + '_tmp.txt'

    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6]

    for cmd in cmds:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode !=0:
            raise Exception(stderr)



def main():

    parser=argparse.ArgumentParser(description="This software extracts inexact HPO terms and generate a gene list.")
    parser.add_argument("-s","--sampleid",help="Sample ID." ,dest="sampleid", type=str, required = True)
    parser.add_argument("-w","--workdir",help="This is the base work directory (folder/directory).",dest="workdir",type=str, required = True)
    parser.add_argument("-j", "--json", help="Medical history JSON file", dest = "json", type=str, required = True)
    parser.add_argument("-m", "--manual", help="Indicate this flag to NOT run CNLP and just use a list of manually curated HPO terms", dest = "manual", action='store_true')
    args=parser.parse_args()


    # run query
    if not args.manual:
        print('[run_CNLP.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Running CNLP on ' + args.sampleid + '...')
        run_clinphen(args)
    else:
        print('[run_CNLP.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Skipping CNLP on ' + args.sampleid + '...')

    print('[run_CNLP.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Gathering additional HPO terms for ' + args.sampleid + '...')
    relatives = queryHPO(args)

    print('[run_CNLP.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Finding clinically relevant genes for ' + args.sampleid + '...')
    getGeneList(args, relatives)




if __name__=="__main__":
    main()




