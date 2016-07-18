#!/bin/bash

module unload ruby
module load ruby

# Deprecated, only needed if you are running RAST.
export RAST_USER=""
export RAST_PASSWORD=""

# Defaults will probably work for these
export PERL5LIB="/usr/bin/perl5.10.1"
export TMP="/sc/orga/scratch/$USER/tmp"
export SMRTANALYSIS="/sc/orga/projects/pacbio/modules/smrtanalysis/2.2.0/install/smrtanalysis_2.3.0.140936"
export SHARED_DIR="/sc/orga/scratch/$USER/shared_dir"
export IGB_DIR="/hpc/users/vanbah01/www/igb"
export REORIENT_FASTA="/sc/orga/projects/InfectiousDisease/reference-db/landmarks/dnaA-reference-sequences.fasta"
export CLUSTER="LSF_PSP"

# If running from interactive1/interactive2, need to run requests through internal HTTP proxy
export HTTP_PROXY="http://proxy.mgmt.hpc.mssm.edu:8123"

# Ensures that the required module files are in MODULEPATH
if [[ ":$MODULEPATH:" != *":/hpc/packages/minerva-mothra/modulefiles:"* ]]; then
    export MODULEPATH="${MODULEPATH:+"$MODULEPATH:"}/hpc/packages/minerva-mothra/modulefiles"
fi
