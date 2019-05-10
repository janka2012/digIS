#!/bin/bash

export INPUT_FASTA
export GENBANK=""
export OUTPUTDIR="digis_output"
export FORMAT="csv"

while getopts f:g:i:o: option
do
case "${option}"
in
f) FORMAT=${OPTARG};;
g) GENBANK=${OPTARG};;
i) INPUT_FASTA=${OPTARG};;
o) OUTPUTDIR=${OPTARG};;
esac
done

DOCKER_DATA=/digis_data
DOCKER_RUN="sudo docker run"
DOCKER_IMAGE="digis"
DOCKER_MOUNT=""
DOCKER_PARAMS=""

[ -z ${INPUT_FASTA} ] && { echo "INPUT_FASTA is unset. This argument is mandatory."; exit 2; } || echo "INPUT_FASTA is set to '$INPUT_FASTA'"
[ -z ${GENBANK} ] && echo "GENBANK is unset" || echo "GENBANK is set to '$GENBANK'"
[ -z ${FORMAT} ] && echo "FORMAT is unset" || echo "FORMAT is set to '$FORMAT'"
[ -z ${OUTPUTDIR} ] && echo "OUTPUTDIR is unset" || echo "OUTPUTDIR is set to '$OUTPUTDIR'"

# mandatory
if [ ! -z ${INPUT_FASTA} ]
then
	input_fasta_fullpath=$(echo $(cd $(dirname "$INPUT_FASTA") && pwd -P)/$(basename "$INPUT_FASTA"))
	DOCKER_MOUNT+=" -v ${input_fasta_fullpath}:${DOCKER_DATA}/$(basename ${INPUT_FASTA})"
	DOCKER_PARAMS+=" -i ${DOCKER_DATA}/${INPUT_FASTA}"
else
	echo "INPUT_FASTA is unset. This argument is mandatory."
	exit 2
fi

# optional
if [ ! -z ${GENBANK} ]
then
	genbank_fullpath=$(echo $(cd $(dirname "$GENBANK") && pwd -P)/$(basename "$GENBANK"))
	DOCKER_MOUNT+=" -v ${genbank_fullpath}:${DOCKER_DATA}/$(basename ${GENBANK})"
	DOCKER_PARAMS+=" -g ${DOCKER_DATA}/${GENBANK}"
fi

DOCKER_PARAMS+=" -f $FORMAT -o $OUTPUTDIR"

INPUT_FASTA=$(basename ${INPUT_FASTA})
GENBANK=$(basename ${GENBANK})

DOCKER_CMD=$DOCKER_RUN$DOCKER_MOUNT" "$DOCKER_IMAGE$DOCKER_PARAMS
echo "$DOCKER_CMD"

sudo docker cp $DOCKER_IMAGE:/digis_output $OUTPUTDIR

