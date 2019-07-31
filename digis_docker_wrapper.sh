#!/bin/bash

###########################################################################################
#
# sudo docker run -v $(pwd)/data/test_data/NC_002608.fasta:/digis_data/NC_002608.fasta
# -v $(pwd)/data/test_data/NC_002608.gb:/digis_data/NC_002608.gb
# -v $(pwd)/digis_output:/digis_output digis
# -i /digis_data/NC_002608.fasta
# -g /digis_data/NC_002608.gb
###########################################################################################

export INPUT_FASTA
export GENBANK
export OUTPUTDIR="digis_output"

while getopts g:i:o: option
do
case "${option}"
in
g) GENBANK=${OPTARG};;
i) INPUT_FASTA=${OPTARG};;
o) OUTPUTDIR=${OPTARG};;
\?) echo "Invalid option: -"$OPTARG"" >&2
    exit 1;;
: ) echo "Option -"$OPTARG" requires an argument." >&2
    exit 1;;
esac
done

DOCKER_DATA="/digis_data"
DOCKER_CONTAINER_NAME="digis_tool"
DOCKER_RUN="sudo docker run --name ${DOCKER_CONTAINER_NAME}"
DOCKER_IMAGE_NAME="digis"
DOCKER_OUTPUTDIR="/digis_output"
DOCKER_MOUNT=""
DOCKER_PARAMS=""

[ -z ${INPUT_FASTA} ] && { echo "INPUT_FASTA is unset. This argument is mandatory."; exit 2; } || echo "INPUT_FASTA is set to '$INPUT_FASTA'"
[ -z ${GENBANK} ] && echo "GENBANK is unset" || echo "GENBANK is set to '$GENBANK'"
[ -z ${OUTPUTDIR} ] && echo "OUTPUTDIR is unset" || echo "OUTPUTDIR is set to '$OUTPUTDIR'"

# mandatory
if [ ! -z ${INPUT_FASTA} ]
then
	input_fasta_fullpath=$(echo $(cd $(dirname "$INPUT_FASTA") && pwd -P)/$(basename "$INPUT_FASTA"))
	DOCKER_MOUNT+=" -v ${input_fasta_fullpath}:${DOCKER_DATA}/$(basename ${INPUT_FASTA})"
	DOCKER_PARAMS+=" -i ${DOCKER_DATA}/$(basename $INPUT_FASTA)"
else
	echo "INPUT_FASTA is unset. This argument is mandatory."
	exit 2
fi

# optional
if [ ! -z ${GENBANK} ]
then
	genbank_fullpath=$(echo $(cd $(dirname "$GENBANK") && pwd -P)/$(basename "$GENBANK"))
	DOCKER_MOUNT+=" -v ${genbank_fullpath}:${DOCKER_DATA}/$(basename ${GENBANK})"
	DOCKER_PARAMS+=" -g ${DOCKER_DATA}/$(basename ${GENBANK})"
fi

# mounting output directory
outputdir_fullpath=$(echo $(cd $(dirname "$OUTPUTDIR") && pwd -P)/$(basename "$OUTPUTDIR"))
DOCKER_MOUNT+=" -v ${outputdir_fullpath}:${DOCKER_OUTPUTDIR}"

# adding script params
DOCKER_PARAMS+=" -o ${DOCKER_OUTPUTDIR}"

# creating full command
DOCKER_CMD=$DOCKER_RUN$DOCKER_MOUNT" "$DOCKER_IMAGE_NAME" python3 digIS_search.py"$DOCKER_PARAMS

# run full command
eval $DOCKER_CMD

# copy output from docker container into host
sudo docker cp ${DOCKER_CONTAINER_NAME}:${DOCKER_OUTPUTDIR} ${outputdir_fullpath}

# remove stopped container
sudo docker container rm $DOCKER_CONTAINER_NAME
