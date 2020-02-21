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
INTERVAR=${5}
BIONANO=${6}
ENZYME=${7}
LINKEDREADSV=${8}

if [[ ! -d ${WORKDIR}/log ]]; then
    mkdir ${WORKDIR}/log
fi


LOGFILE=${WORKDIR}/log/`date '+%Y-%m-%d'`_pipeline_run_${SAMPLEID}.log

if [[ -f ${LOGFILE} ]]; then
    rm ${LOGFILE}
fi


## main pipeline
pipeline(){
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 


if [[ ! -d ${WORKDIR}/results ]];
    then
    mkdir ${WORKDIR}/results
fi

if [[ ! -d ${WORKDIR}/results/${SAMPLEID} ]];
    then
    mkdir ${WORKDIR}/results/${SAMPLEID}
fi

additional_var=''
if [[ ${BIONANO} = true ]]; then
    additional_var+=" -b -e ${ENZYME}"
fi
if [[ ${LINKEDREADSV} = true ]]; then
    additional_var+=" -l"
fi

# extract exact and inexact HPO terms then run CNLP
python3.6 ${WORKDIR}/scripts/run_CNLP.py -s ${SAMPLEID} -w ${WORKDIR} -d ${DATABASE} -j ${JSON}

# identify variants and rank them 
python3.6 ${WORKDIR}/scripts/run_clinical_interpretor.py -s ${SAMPLEID} -w ${WORKDIR} -d ${DATABASE} -i ${INTERVAR} ${additional_var}

} # end of pipeline

pipeline 2>&1 | tee $LOGFILE
