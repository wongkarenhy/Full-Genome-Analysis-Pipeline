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
OUT=${2}
SAMPLEID=${3}

# extract from notes/medical history
## features/notes
cat $JSON | jq-linux64 '.features[].notes' | grep -v null | sed -e 's/$/\./g' > $OUT/$SAMPLEID_tmp.txt
cat $JSON | jq-linux64 '.notes.medical_history' >> $OUT/$SAMPLEID_tmp.txt
cat $JSON | jq-linux64 '.features[].id' | tr -d '"' > ${OUT}/${SAMPLEID}_hpo_manual.txt

# clean up the file
cat $OUT/$SAMPLEID_tmp.txt | tr -d '[]"\n' | tr '\' '!' | sed 's/!r!n/\./g' > $OUT/$SAMPLEID.txt

# extract HPO terms
clinphen $OUT/$SAMPLEID.txt | awk -F'\t' '{print $1}' | tail -n +2 > ${OUT}/${SAMPLEID}_hpo_exact.txt

# extract inexact HPO terms
python3 extract_inexact_terms.py -s $SAMPLEID 

# remove intermediate files
rm $OUT/$SAMPLEID_tmp.txt
