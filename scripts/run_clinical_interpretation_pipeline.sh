#!/bin/bash

# exit immediately upon error 
set -e

while getopts 'j:w:s:i:b:e:l:t:f:m:r:a:' OPTION; do
  case "$OPTION" in
    j)
      JSON="$OPTARG"
      ;;
    w)
      WORKDIR="$OPTARG"
      ;;
    s)
      SAMPLEID="$OPTARG"
      ;;
    i)
      INTERVAR="$OPTARG"
      ;;
    b)
      BIONANO="$OPTARG"
      [[ ! $BIONANO =~ true|false ]] && {
          echo "Bionano flag can be either true or false"
          exit 1
      }

      ;;
    e)
      ENZYME="$OPTARG"
      [[ ! $ENZYME =~ BspQI|DLE|None ]] && {
          echo "Incorrect options provided in ENZYME"
          exit 1
      }

      ;;
    l)
      LINKEDREADSV="$OPTARG"
      [[ ! $LINKEDREADSV =~ true|false ]] && {
          echo "LINKEDREADSV can be either true or false"
          exit 1
      }

      ;;
    t)
      TYPE="$OPTARG"
      [[ ! $TYPE =~ trio|singleton ]] && {
	  echo "Incorrect options provided in TYPE"
          exit 1
      }
      ;;
    f)
      FATHERVCF="$OPTARG"
      ;;
    m)
      MOTHERVCF="$OPTARG"
      ;;
    r)
      REF="$OPTARG"
      [[ ! $REF =~ hg19|hg38 ]] && {
          echo "Incorrect options provided in REF"
          exit 1
      }

      ;;
    a)
      ARTIFACT="$OPTARG"
      ;;

    ?)
      echo "script usage: $(basename $0) [-j path_to_json/None] [-w work_dir] [-s sample_id] [-i path_to_intervar] [-b true/false] [-e DLE/BspQI/None] [-l true/false] [-t trio/singleton] [-f father_SNP_vcf_file_path or None if singleton] [-m mother_SNP_vcf_file_path or None if singleton] [-r hg19/hg38] [-a path_to_custom_artifact_file or None]" >&2
      exit 1
      ;;
  esac
  
done
shift $(( OPTIND - 1 ))

if [[ ! -d ${WORKDIR}/log ]]; then
    mkdir ${WORKDIR}/log
fi


LOGFILE=${WORKDIR}/log/`date '+%Y-%m-%d'`_pipeline_run_${SAMPLEID}.log

if [[ -f ${LOGFILE} ]]; then
    rm ${LOGFILE}
fi


# main pipeline
pipeline(){
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START:  ${0} -j $JSON -w $WORKDIR -s $SAMPLEID -i $INTERVAR -b $BIONANO -e $ENZYME -l $LINKEDREADSV -t $TYPE -f $FATHERVCF -m $MOTHERVCF -r $REF -a $ARTIFACT"

if [[ ! -d ${WORKDIR}/results ]];
    then
    mkdir ${WORKDIR}/results
fi

if [[ ! -d ${WORKDIR}/results/${SAMPLEID} ]];
    then
    mkdir ${WORKDIR}/results/${SAMPLEID}
fi

if [[ ! -d ${WORKDIR}/results/${SAMPLEID}/confident_set ]]; then
    mkdir ${WORKDIR}/results/${SAMPLEID}/confident_set
fi
if [[ ! -d ${WORKDIR}/results/${SAMPLEID}/misc ]]; then
    mkdir ${WORKDIR}/results/${SAMPLEID}/misc
fi

additional_var=''
if [[ ${BIONANO} = true ]]; then
    additional_var+=" -b -e ${ENZYME}"
fi
if [[ ${LINKEDREADSV} = true ]]; then
    additional_var+=" -l"
fi
if [[ ${TYPE} == 'singleton' ]]; then
    additional_var+=" -S"
else
    additional_var+=" -f ${FATHERVCF} -m ${MOTHERVCF}"
fi

additional_var_CNLP=''
if [[ ${JSON} = None ]]; then
    additional_var_CNLP+=" -m"
fi


# extract exact and inexact HPO terms then run CNLP
python3.6 ${WORKDIR}/scripts/run_CNLP.py -s ${SAMPLEID} -w ${WORKDIR} -j ${JSON} ${additional_var_CNLP}

# identify variants and rank them 
python3.6 ${WORKDIR}/scripts/run_clinical_interpretor.py -s ${SAMPLEID} -w ${WORKDIR} -i ${INTERVAR} -r ${REF} -a ${ARTIFACT} ${additional_var}

# generate an html report
python3.6 ${WORKDIR}/scripts/generate_report.py -s ${SAMPLEID} -w ${WORKDIR}/results/${SAMPLEID}/confident_set/

} # end of pipeline
pipeline 2>&1 | tee $LOGFILE

