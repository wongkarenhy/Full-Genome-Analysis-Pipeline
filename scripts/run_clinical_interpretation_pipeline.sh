#!/bin/bash

# exit immediately upon error 
set -e

while getopts 'j:w:s:i:b:e:l:t:f:m:r:a:x:M:' OPTION; do
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
          echo "ERROR! Pipeline aborted! Bionano flag has to be either true or false"
          exit 1
      }

      ;;
    e)
      ENZYME="$OPTARG"
      [[ ! $ENZYME =~ BspQI|DLE|None ]] && {
          echo "ERROR! Pipeline aborted! Incorrect option provided in ENZYME"
          exit 1
      }

      ;;
    l)
      LINKEDREADSV="$OPTARG"
      [[ ! $LINKEDREADSV =~ true|false ]] && {
          echo "ERROR! Pipeline aborted! LINKEDREADSV has to be either true or false"
          exit 1
      }

      ;;
    t)
      TYPE="$OPTARG"
      [[ ! $TYPE =~ trio|singleton|duo ]] && {
	  echo "ERROR! Pipeline aborted! Incorrect option provided in TYPE"
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
          echo "ERROR! Pipeline aborted! Incorrect option provided in REF"
          exit 1
      }

      ;;
    x)
      XLINK="$OPTARG"
      [[ ! $XLINK =~ true|false ]] && {
          echo "ERROR! Pipeline aborted! XLINK has to be either true or false"
          exit 1
      }

      ;;
    a)
      ARTIFACT="$OPTARG"
      ;;
    M)
      MAF="$OPTARG"
      ;;

    ?)
      echo "script usage: $(basename $0) [-j path_to_json/None] [-w work_dir] [-s sample_id] [-i path_to_intervar] [-b true/false] [-e DLE/BspQI/None] [-l true/false] [-t trio/singleton] [-f father_SNP_vcf_file_path or None if singleton] [-m mother_SNP_vcf_file_path or None if singleton] [-r hg19/hg38] [-a path_to_custom_artifact_file or None] [-x true/false] [-M 0-1]" >&2
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
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START:  ${0} -j $JSON -w $WORKDIR -s $SAMPLEID -i $INTERVAR -b $BIONANO -e $ENZYME -l $LINKEDREADSV -t $TYPE -f $FATHERVCF -m $MOTHERVCF -r $REF -a $ARTIFACT -x $XLINK -M $MAF" 

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

# If type is trio, both FATHERVCF and MOTHERVCF cannot be None
# If type is duo, either FATHERVCF or MOTHERVCF must be specified 
if [[ ${TYPE} == 'trio' && (${FATHERVCF} == 'None' || ${MOTHERVCF} == 'None') ]]; then
    echo 'ERROR! Pipeline aborted! If type is trio, father and mother VCF files must both be specified.'
    exit 1
elif [[ ${TYPE} == 'duo' && ${FATHERVCF} == 'None' && ${MOTHERVCF} == 'None' ]]; then
    echo 'ERROR! Pipeline aborted! If type is duo, either father or mother VCF file must be specified.'
    exit 1
elif [[ ${TYPE} == 'duo' && ${FATHERVCF} != 'None' && ${MOTHERVCF} != 'None' ]]; then
    echo 'ERROR! Pipeline aborted! If type is duo, either father or mother VCF file has to be None.'
    exit 1
elif [[ ${TYPE} == 'singleton' && (${FATHERVCF} != 'None' || ${MOTHERVCF} != 'None') ]]; then
    echo 'ERROR! Pipeline aborted! If type is singleton, both father or mother VCF files have to be None.'
    exit 1
fi


additional_var=''
if [[ ${BIONANO} == 'true' ]]; then
    additional_var+=" -b -e ${ENZYME}"
fi
if [[ ${LINKEDREADSV} == 'true' ]]; then
    additional_var+=" -l"
fi

# Parse type if duo
if [[ ${TYPE} == 'duo' ]]; then

    if [[ ${FATHERVCF} != 'None' ]]; then
	additional_var+=" -F"
    else
	additional_var+=" -M"
    fi
fi

# Check if the pipeline should generate xlink small variants file
if [[ ${XLINK} == 'true' ]]; then
    additional_var+=" -X"
fi

additional_var_CNLP=''
if [[ ${JSON} == 'None' ]]; then
    additional_var_CNLP+=" -m"
fi

# Check if the intervar output file is present and not empty
if [[ -s ${INTERVAR}/example/${SAMPLEID}.hg38_multianno.txt.intervar ]]; then
    grep -w -E 'frameshift|nonframeshift|nonsynonymous|stopgain|stoploss|splicing' ${INTERVAR}/example/${SAMPLEID}.${REF}_multianno.txt.intervar | awk -F'\t' -v maf=$MAF '$15<=maf' > ${INTERVAR}/example/${SAMPLEID}.${REF}_multianno.txt.intervar.FINAL
    awk -F'\t' '{print $6}' ${INTERVAR}/example/${SAMPLEID}.${REF}_multianno.txt.intervar.FINAL | sort -u > ${INTERVAR}/example/${SAMPLEID}_smallVariant_geneList.txt
else
    echo "${INTERVAR}/example/${SAMPLEID}.hg38_multianno.txt.intervar is either missing or empty"
    exit 1
fi

# extract exact and inexact HPO terms then run CNLP
python3.6 ${WORKDIR}/scripts/run_CNLP.py -s ${SAMPLEID} -w ${WORKDIR} -j ${JSON} ${additional_var_CNLP}

# identify variants and rank them
python3.6 ${WORKDIR}/scripts/run_clinical_interpretor.py -s ${SAMPLEID} -w ${WORKDIR} -i ${INTERVAR} -r ${REF} -a ${ARTIFACT} -t ${TYPE} -f ${FATHERVCF} -m ${MOTHERVCF} ${additional_var}

# generate an html report
python3.6 ${WORKDIR}/scripts/generate_report.py -s ${SAMPLEID} -w ${WORKDIR}/results/${SAMPLEID}/confident_set/

} # end of pipeline
pipeline 2>&1 | tee $LOGFILE

