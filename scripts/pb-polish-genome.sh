#!/bin/sh

# 01.10.2019 12:46:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# Process command line arguments
SUBREADS=""
REFERENCE=""
OUTPREFIX="out"
ALGORITHM="arrow"
CPUS=8
while getopts "i:r:o:c:" opt; do
   case $opt in
     
   i)
      SUBREADS="$OPTARG"
      ;;
   
   r)
      REFERENCE="$OPTARG"
      ;;
   
   o)
      OUTPREFIX="$OPTARG"
      ;;

   c) 
      CPUS="$OPTARG"
      [[ ! $CPUS =~ ^[0-9]+$ ]] && {
         echo "Incorrect CPU option provided"
         exit 1
      }
      ;;
   
   a)
      ALGORITHM="$OPTARG"
      ;;
   
   *)
      echo "Incorrect options provided"
      exit 1
      ;;

   esac
done

# Check arguments and produce a help message
if [ -z "$SUBREADS" ];
then
  cat << EOF

   Usage: pb-polish-genome.sh -i <subread.bam> -r <reference-fasta> -o <output prefix>

   Arguments:
    -i <string>
      Subread bam file. Must be created by sequel2 pipeline or converted
      from bax.h5 using the pbconda bax2bam utility.
    -p <string>
      Fasta file with reference genome
    -o <string>
      Output prefix
    -c <integer>
      Number of CPUs. Default: 8
    -help
      This help message

EOF
  exit 0
fi

# Test case
#SUBREADS="ER15890_3A_028012.subreads.bam"
#REFERENCE="S_aureus_ER15890_3A_028012.genome.fasta"
#OUTPREFIX="polish_ER15890_3A_028012"
#CPUS=8

######################
# SET UP ENVIRONMENT #
######################

# Load anaconda module
module purge all
unset PYTHONPATH
unset PERL5LIB
unset R_LIBS
module load anaconda2
module load zlib

# Create conda environment if it doesn't already exist
checkenv=`conda info --envs | grep pbpolish | wc -l`
if [ "checkenv" == "0" ]
then
   conda create -n pbpolish pbmm2 pbcore pbcoretools pbgcpp bax2bam genomicconsensus samtools
else
   conda deactivate
   source activate pbpolish
fi

###########
# ROUND 1 #
###########

echo "Subread file    : ${SUBREADS}"
echo "Reference fasta : ${REFERENCE}"
echo "Output prefix   : ${OUTPREFIX}"
echo "Algorithm       : ${ALGORITHM}"
echo "CPUs            : ${CPUS}"

# Align reads to reference
pbmm2 align --sort -j ${CPUS} -J 2 ${REFERENCE} ${SUBREADS} ${OUTPREFIX}_align0.bam 

# Use either gcpp (for new sequel chemistry) or the arrow algorithm
if [ "$ALGORITHM" == "arrow" ]
then
   pbindex ${OUTPREFIX}_align0.bam
   samtools faidx ${REFERENCE}
   variantCaller -j ${CPUS} --algorithm arrow -r ${REFERENCE} -o ${OUTPREFIX}_polish0.fasta -o ${OUTPREFIX}_polish0.gff -o ${OUTPREFIX}_polish0.vcf ${OUTPREFIX}_align0.bam
else
   gcpp -j ${CPUS} -r ${REFERENCE} -o ${OUTPREFIX}_polish0.fasta ${OUTPREFIX}_align0.bam
fi

###########
# ROUND 2 #
###########

# Align reads to reference
pbmm2 align --sort -j ${CPUS} -J 2 ${OUTPREFIX}_polish0.fasta ${SUBREADS} ${OUTPREFIX}_align1.bam 

# Use either gcpp (for new sequel chemistry) or the arrow algorithm
if [ "$ALGORITHM" == "arrow" ]
then
   pbindex ${OUTPREFIX}_align1.bam 
   samtools faidx ${OUTPREFIX}_polish0.fasta
   variantCaller -j ${CPUS} --algorithm arrow -r ${OUTPREFIX}_polish0.fasta -o ${OUTPREFIX}_polish1.fasta -o ${OUTPREFIX}_polish1.gff -o ${OUTPREFIX}_polish1.vcf ${OUTPREFIX}_align1.bam
else
   gcpp -j ${CPUS} -r ${OUTPREFIX}_polish0.fasta -o ${OUTPREFIX}_polish1.fasta ${OUTPREFIX}_align1.bam
fi
