#!/bin/sh

# 17.03.2019 12:56:56 EDT
# Harm van Bakel <hvbakel@gmail.com>

module purge all
module load python/2.7.6
module load py_packages/2.7
module load blat
module load bioperl
module load ucsc-utils/2015-04-07
module load openssl/1.0.2

for i in ER*prokka
do
   cwd=`pwd`
   name=`basename $i _prokka`
   cd $i
   rast2igb.pl -f ${i}.gbk -g E_faecalis_${name} -i /sc/orga/projects/InfectiousDisease/igb -r ~/opt/pathogendb-pipeline/
   igb2pathogendb.pl -i /sc/orga/projects/InfectiousDisease/igb/E_faecalis_${name}
   cd $cwd
done
