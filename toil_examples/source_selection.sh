#!/bin/bash

#NOTE: works only with TOIL>6.0.0

#### UPDATE THESE ####

export TOIL_SLURM_ARGS="--export=ALL --job-name phasediff -p normal -t 4:00:00"
SING_BIND="/project,/project/lofarvwf/Software,/project/lofarvwf/Share,/project/lofarvwf/Public,/home/lofarvwf-jdejong"
VLBI_SCRIPTS="/project/lofarvwf/Software/vlbi/scripts/"

CWL_WORKFLOW=/project/lofarvwf/Software/lofar_helpers/source_selection/pre_selection/ddcal_pre_selection.cwl
YAML=input.yaml
VENV=/home/lofarvwf-jdejong/venv

######################

# set up singularity
SIMG=vlbi-cwl.sif
mkdir -p singularity
wget https://lofar-webdav.grid.sara.nl/software/shub_mirror/tikk3r/lofar-grid-hpccloud/amd/flocs_v4.5.0_znver2_znver2_aocl4_cuda.sif -O singularity/$SIMG
mkdir -p singularity/pull
cp singularity/$SIMG singularity/pull/$SIMG

CONTAINERSTR=$(singularity --version)
if [[ "$CONTAINERSTR" == *"apptainer"* ]]; then
  export APPTAINER_CACHEDIR=$PWD/singularity
  export APPTAINER_TMPDIR=$APPTAINER_CACHEDIR/tmp
  export APPTAINER_PULLDIR=$APPTAINER_CACHEDIR/pull
  export APPTAINER_BIND=$SING_BIND
  export APPTAINERENV_PYTHONPATH='${VLBI_SCRIPTS}:$PYTHONPATH'
else
  export SINGULARITY_CACHEDIR=$PWD/singularity
  export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
  export SINGULARITY_PULLDIR=$SINGULARITY_CACHEDIR/pull
  export SINGULARITY_BIND=$SING_BIND
  export SINGULARITYENV_PYTHONPATH='${VLBI_SCRIPTS}:$PYTHONPATH'
fi

export CWL_SINGULARITY_CACHE=$APPTAINER_CACHEDIR
export TOIL_CHECK_ENV=True

# make folder for running toil
WORKDIR=$PWD/workdir
OUTPUT=$PWD/outdir
JOBSTORE=$PWD/jobstore
LOGDIR=$PWD/logs
TMPD=$PWD/tmpdir

mkdir -p ${TMPD}_interm
mkdir -p $WORKDIR
mkdir -p $OUTPUT
mkdir -p $LOGDIR

source ${VENV}/bin/activate

# run toil
toil-cwl-runner \
--no-read-only \
--retryCount 0 \
--singularity \
--disableCaching \
--writeLogsFromAllJobs True \
--logFile full_log.log \
--writeLogs ${LOGDIR} \
--outdir ${OUTPUT} \
--tmp-outdir-prefix ${TMPD}/ \
--jobStore ${JOBSTORE} \
--workDir ${WORKDIR} \
--coordinationDir ${OUTPUT} \
--tmpdir-prefix ${TMPD}_interm/ \
--disableAutoDeployment True \
--bypass-file-store \
--preserve-entire-environment \
--batchSystem slurm \
${CWL_WORKFLOW} ${YAML}

#--cleanWorkDir never \ --> for testing

deactivate
