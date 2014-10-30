#!/bin/bash

module unload ruby
# module unload python

module load ruby
# module load python/2.7.3
# module load py_packages
# module unload gcc
# module load smrtpipe/2.2.0

export PERL5LIB="/usr/bin/perl5.10.1"
export TMP="/sc/orga/scratch/$USER"
export SMRTPIPE="/sc/orga/projects/InfectiousDisease/smrtpipe"
export SMRTANALYSIS="/hpc/packages/minerva-mothra/smrtanalysis/2.2.0/ROOT"
