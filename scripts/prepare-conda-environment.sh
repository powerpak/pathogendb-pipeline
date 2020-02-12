#!/bin/sh

# 11.02.2020 19:28:53 EST
# Harm van Bakel <hvbakel@gmail.com>

# Set anaconda environment
module purge all
unset PYTHONPATH
unset PERL5LIB
unset R_LIBS
module load anaconda2
module load zlib

# Prepare prokka conda environment
conda create -n prokka -c bioconda prokka=1.14.5
wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz -O linux64.tbl2asn.gz 
gunzip linux64.tbl2asn.gz
mv linux64.tbl2asn ~/.conda/envs/prokka/bin/tbl2asn
chmod +x ~/.conda/envs/prokka/bin/tbl2asn

# Prepare pbpolish conda environment
conda create -n pbpolish pbmm2 pbcore pbcoretools pbgcpp bax2bam genomicconsensus samtools
