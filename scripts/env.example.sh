#!/bin/bash

module unload ruby
module load ruby

# Must be configured by the end user
export RAST_USER=""
export RAST_PASSWORD=""

# Defaults will probably work for these
export PERL5LIB="/usr/bin/perl5.10.1"
export TMP="/sc/orga/scratch/$USER/tmp"
export SMRTPIPE="/sc/orga/projects/InfectiousDisease/smrtpipe"
export SMRTANALYSIS="/sc/orga/projects/pacbio/modules/smrtanalysis/2.2.0/install/smrtanalysis_2.3.0.140936"
export SHARED_DIR="/sc/orga/scratch/$USER/shared_dir"

# Ensures that the required module files are in MODULEPATH
if [[ ":$MODULEPATH:" != *":/hpc/packages/minerva-mothra/modulefiles:"* ]]; then
    export MODULEPATH="${MODULEPATH:+"$MODULEPATH:"}/hpc/packages/minerva-mothra/modulefiles"
fi
