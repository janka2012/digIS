#!/bin/bash

export INPUT_FASTA
export GENBANK=""
export OUTPUTDIR="digis_output"
export FORMAT="csv"

usage()
{
  echo "Usage: $0 -i INPUT_FASTA [ -g GENBANK ] [-f FORMAT] [-o OUTPURDIR]"
}

show_help() {
usage
cat << EOF
-h	display this help and exit.
-f	specify output format, csv or gff, default=csv.
-g	(multi)genbank annotation file, optional.
-i	(multi)fasta file, required.
-o	output directory name, default=digis_output.
EOF
}

while getopts :hf:g:i:o: option
do
case "${option}"
in
h | \? | * )
	show_help
	exit 1
	;;
f )
	FORMAT=${OPTARG}
	;;
g )
	GENBANK=${OPTARG}
	;;
i )
	INPUT_FASTA=${OPTARG}
	;;
o )
	OUTPUTDIR=${OPTARG}
	;;
: )
	echo "Invalid option: $OPTARG requires an argument" 1>&2
	;;
esac
done

DOCKER_DATA="/digis_data"
UNIQ_ID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 10 | head -n 1)
DOCKER_CONTAINER_NAME="digis_${UNIQ_ID}"
DOCKER_RUN="sudo docker run --name ${DOCKER_CONTAINER_NAME} "
DOCKER_IMAGE_NAME="digis"
DOCKER_OUTPUTDIR="/digis_output"
DOCKER_MOUNT=""
DOCKER_PARAMS=""

# mandatory
if [ ! -z ${INPUT_FASTA} ]
then
	echo "INPUT_FASTA is set to ${INPUT_FASTA}"
	input_fasta_fullpath=$(echo $(cd $(dirname ${INPUT_FASTA}) && pwd -P)/$(basename ${INPUT_FASTA}))
	DOCKER_MOUNT="${DOCKER_MOUNT} -v ${input_fasta_fullpath}:${DOCKER_DATA}/$(basename ${INPUT_FASTA})"
	DOCKER_PARAMS="${DOCKER_PARAMS} -i ${DOCKER_DATA}/$(basename ${INPUT_FASTA})"

else
	echo "INPUT_FASTA is unset. This argument is mandatory." >&2
	exit 1
fi

# optional
if [ ! -z ${GENBANK} ]
then
	echo "GENBANK is set to ${GENBANK}"
	genbank_fullpath=$(echo $(cd $(dirname "$GENBANK") && pwd -P)/$(basename "$GENBANK"))
	DOCKER_MOUNT="${DOCKER_MOUNT} -v ${genbank_fullpath}:${DOCKER_DATA}/$(basename ${GENBANK})"
	DOCKER_PARAMS="${DOCKER_PARAMS} -g ${DOCKER_DATA}/$(basename $GENBANK)"
else
	echo "No GenBank annotation file provided."
fi

[ -z ${FORMAT} ] && echo "FORMAT is unset" || echo "FORMAT is set to '$FORMAT'"
[ -z ${OUTPUTDIR} ] && echo "OUTPUTDIR is unset" || echo "OUTPUTDIR is set to '$OUTPUTDIR'"

# mounting output directory
outputdir_fullpath=$(echo $(cd $(dirname "$OUTPUTDIR") && pwd -P)/$(basename "$OUTPUTDIR"))
DOCKER_MOUNT="${DOCKER_MOUNT} -v ${outputdir_fullpath}:${DOCKER_OUTPUTDIR}"

# adding script params
DOCKER_PARAMS="${DOCKER_PARAMS} -f $FORMAT -o ${DOCKER_OUTPUTDIR}"

# creating full command
DOCKER_CMD=${DOCKER_RUN}${DOCKER_MOUNT}" "${DOCKER_IMAGE_NAME}" python3 digIS_search.py"${DOCKER_PARAMS}

# run full command
eval $DOCKER_CMD

# remove stopped container
sudo docker container rm $DOCKER_CONTAINER_NAME
