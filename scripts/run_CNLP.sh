#!/bin/bash                                                                                                                                                                                                 

# exit immediately upon error                                                                                                                                                                               
set -e

usage() {
    NAME=$(basename $0)
    cat <<EOF
Usage:                                                                                                                                                                                                      
  ${NAME} [path_to_json] [out_path] [sample_id]                                                                                                                                                             
This software extracts HPO terms from json files, searches for inexact terms, and outputs the associated genes.

EOF
}

# declare variable                                                                                                                                                                                          
JSON=${1}
WORKDIR=${2}
SAMPLEID=${3}
DATABASE=${4}

if [[ ! -d ${WORKDIR}/results ]];
    then
    mkdir ${WORKDIR}/results
fi

if [[ ! -d ${WORKDIR}/results/${SAMPLEID} ]];
    then
    mkdir ${WORKDIR}/results/${SAMPLEID}
fi

# extract from notes/medical history                                                                                                                                                                        
## features/notes                                                                                                                                                                                           
cat $JSON | jq-linux64 '.features[].notes' | grep -v null | sed -e 's/$/\./g' > ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}_tmp.txt
cat $JSON | jq-linux64 '.notes.medical_history' >> ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}_tmp.txt
cat $JSON | jq-linux64 '.features[].id' | tr -d '"' > ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}_hpo_manual.txt

# clean up the file                                                                                                                                                                                         
cat ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}_tmp.txt | tr -d '[]"\n' | tr '\' '!' | sed 's/!r!n/\./g' > ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}.txt

# extract HPO terms                                                                                                                                                                                         
clinphen ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}.txt | awk -F'\t' '{print $1}' | tail -n +2 > ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}_hpo_exact.txt

# extract inexact HPO terms                                                                                                                                                                                 
python3.6 ${WORKDIR}/scripts/extract_inexact_terms.py -s ${SAMPLEID} -w ${WORKDIR} -d ${DATABASE}

# remove intermediate files                                                                                                                                                                                 
rm ${WORKDIR}/results/${SAMPLEID}/${SAMPLEID}_tmp.txt

