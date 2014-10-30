#!/bin/bash

module load ruby
module load python/2.7.6
module load py_packages

export PERL5LIB="/usr/bin/perl5.10.1"
export TMP="/sc/orga/scratch/$USER"
export SMRTPIPE="/sc/orga/projects/InfectiousDisease/smrtpipe"
export SMRTANALYSIS="/hpc/packages/minerva-mothra/smrtanalysis/2.2.0/ROOT"
