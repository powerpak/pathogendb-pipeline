#!/usr/bin/env python

import sys
import subprocess
import os
import argparse

def filter_vcf(in_file, out_file):
    with open(in_file) as vcf, open(out_file, 'w') as out:
        for line in vcf:
            if line.startswith('#'):
                out.write(line)
            else:
                chrom, pos, id, ref, alt, qual, filt, info, form, unk = line.split()
                #if len(ref) != len(alt):
                for i in info.split(';'):
                    if i.startswith('AF1='):
                        if float(i.split('=')[1]) > 0.9:
                            out.write(line)
                            continue


def correct_regions(fasta_file, read_file, coverage_file, working_dir, out_file, read_file_2, read_length=100):
    low_cov = {}
    with open(coverage_file) as cov:
        cov_dict = {}
        cov_vals = []
        for line in cov:
            ref, pos, cov = line.split()
            if not ref in cov_dict:
                cov_dict[ref] = []
            cov_dict[ref].append(int(cov))
            cov_vals.append(int(cov))
        cov_vals.sort()
        median_cov = cov_vals[len(cov_vals)/2]
        cov_cutoff = median_cov / 8
        sys.stdout.write('Using a coverage cutoff of ' + str(cov_cutoff))
        for i in cov_dict: # find regions in assembly wiht low coverage
            low_cov[i] = []
            minval = None
            maxval = None
            for pos, j in enumerate(cov_dict[i]):
                if j < cov_cutoff:
                    if minval is None:
                        minval = pos - read_length
                    maxval = pos
                else:
                    if not maxval is None and pos > maxval + 2 * read_length:
                        low_cov[i].append((minval, maxval + read_length))
                        minval = None
                        maxval = None
    with open(fasta_file) as fasta:
        seqDict = {}
        for line in fasta:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                seqDict[name] = ''
            else:
                seqDict[name] += line.rstrip()
    split_seq = {}

    for i in low_cov:
        last_pos = 0
        split_seq[i] = []
        for j in low_cov[i]:
            split_seq[i].append(seqDict[i][last_pos:j[0]])
            split_seq[i].append(seqDict[i][j[0]:j[1]])
            last_pos = j[1]
        split_seq[i].append(seqDict[i][last_pos:])
    with open(working_dir + '/ref.fa', 'w') as ref:
        for i in split_seq:
            for num, j in enumerate(split_seq[i]):
                if num % 2 == 1:
                    ref.write('>' + i + '_' + str(num) + '\n')
                    for k in range(0, len(j), 60):
                        ref.write(j[k:k+60] + '\n')
    subprocess.Popen('bwa index ' + working_dir + '/ref.fa', shell=True).wait()
    if read_file_2 is None:
        subprocess.Popen('bwa mem -t 4 ' + working_dir + '/ref.fa ' + read_file + ' > ' + working_dir + '/ref.aln.sam', shell=True).wait()
    else:
        subprocess.Popen('bwa mem -t 4 ' + working_dir + '/ref.fa ' + read_file + ' ' + read_file_2 + ' > ' + working_dir + '/ref.aln.sam', shell=True).wait()
    subprocess.Popen('samtools faidx ' + working_dir + '/ref.fa ', shell=True).wait()
    subprocess.Popen('samtools view -bS ' + working_dir + '/ref.aln.sam > ' + working_dir + '/ref.aln.bam', shell=True).wait()
    subprocess.Popen('samtools sort ' + working_dir + '/ref.aln.bam ' + working_dir + '/ref.sort', shell=True).wait()
    subprocess.Popen('samtools index ' + working_dir + '/ref.sort.bam', shell=True).wait()
    subprocess.Popen('samtools mpileup -L100000 -d100000 -uf "' + working_dir + '/ref.fa" ' + working_dir + \
                    '/ref.sort.bam | bcftools call -cv -Ob > "' + working_dir + '/ref.bcf"', shell=True).wait()
    subprocess.Popen('bcftools view "' + working_dir + '/ref.bcf" | vcfutils.pl varFilter -w 0 -W 0 > "' + working_dir + '/ref.vcf"', shell=True).wait()
    filter_vcf(working_dir + '/ref.vcf', working_dir + '/ref2.vcf')
    subprocess.Popen('bgzip -c "' + working_dir + '/ref2.vcf" > "' + working_dir + '/ref2.vcf.gz"', shell=True).wait()
    subprocess.Popen('tabix -p vcf "' + working_dir + '/ref2.vcf.gz"', shell=True).wait()
    subprocess.Popen('cat "' + working_dir + '/ref.fa" | vcf-consensus "' + working_dir + '/ref2.vcf.gz" > "' + working_dir + '/new_ref.fasta"', shell=True).wait()
    with open(working_dir + '/new_ref.fasta') as fasta:
        for line in fasta:
            if line.startswith('>'):
                name = '_'.join(line.rstrip()[1:].split('_')[:-1])
                num = int(line.rstrip().split('_')[-1])
                split_seq[name][num] = ''
            else:
                split_seq[name][num] += line.rstrip()
    subprocess.Popen('genomeCoverageBed -d -ibam ' + working_dir + '/ref.sort.bam -g ' + working_dir + '/new_ref.fasta > ' + working_dir + '/ref.cov', shell=True).wait()
    redo = set()
    with open(working_dir + '/ref.cov') as cov:
        cov_dict = {}
        for line in cov:
            ref, pos, cov = line.split()
            if not ref in cov_dict:
                cov_dict[ref] = []
            cov_dict[ref].append(int(cov))
        for i in cov_dict:
            for j in cov_dict[i][read_length-10:-read_length+10]:
                if j < cov_cutoff:
                    redo.add(i)
                    continue
    for i in redo:
        with open(working_dir + '/ref.fa', 'w') as ref:
            name = '_'.join(i.split('_')[:-1])
            num = int(i.split('_')[-1])
            ref.write('>' + i + '\n')
            ref.write(split_seq[name][num])
        subprocess.Popen('bwa index ' + working_dir + '/ref.fa', shell=True).wait()
        if read_file_2 is None:
            subprocess.Popen('bwa mem -t 4 ' + working_dir + '/ref.fa ' + read_file + ' > ' + working_dir + '/ref.aln.sam', shell=True).wait()
        else:
            subprocess.Popen('bwa mem -t 4 ' + working_dir + '/ref.fa ' + read_file + ' ' + read_file_2 + ' > ' + working_dir + '/ref.aln.sam', shell=True).wait()
        subprocess.Popen('samtools faidx ' + working_dir + '/ref.fa ', shell=True).wait()
        subprocess.Popen('samtools view -bS ' + working_dir + '/ref.aln.sam > ' + working_dir + '/ref.aln.bam', shell=True).wait()
        subprocess.Popen('samtools sort ' + working_dir + '/ref.aln.bam ' + working_dir + '/ref.sort', shell=True).wait()
        subprocess.Popen('samtools index ' + working_dir + '/ref.sort.bam', shell=True).wait()
        subprocess.Popen('samtools mpileup -L100000 -d100000 -uf "' + working_dir + '/ref.fa" ' + working_dir + \
                        '/ref.sort.bam | bcftools call -cv -Ob > "' + working_dir + '/ref.bcf"', shell=True).wait()
        subprocess.Popen('bcftools view "' + working_dir + '/ref.bcf" | vcfutils.pl varFilter -w 0 -W 0 > "' + working_dir + '/ref.vcf"', shell=True).wait()
        filter_vcf(working_dir + '/ref.vcf', working_dir + '/ref2.vcf')
        subprocess.Popen('bgzip -c "' + working_dir + '/ref2.vcf" > "' + working_dir + '/ref2.vcf.gz"', shell=True).wait()
        subprocess.Popen('tabix -p vcf "' + working_dir + '/ref2.vcf.gz"', shell=True).wait()
        subprocess.Popen('cat "' + working_dir + '/ref.fa" | vcf-consensus "' + working_dir + '/ref2.vcf.gz" > "' + working_dir + '/new_ref.fasta"', shell=True).wait()
        with open(working_dir + '/new_ref.fasta') as fasta:
            seq = ''
            for line in fasta:
                if not line.startswith('>'):
                    seq += line.rstrip()
        split_seq[name][num] = seq
    with open(out_file, 'w') as out:
        for i in split_seq:
            out.write('>' + i + '\n')
            seq = ''.join(split_seq[i])
            for j in range(0, len(seq), 80):
                out.write(seq[j:j+80] + '\n')



parser = argparse.ArgumentParser(prog='Fix_repeats_ill.py', formatter_class=argparse.RawDescriptionHelpFormatter, description='''
fix_repeats_ill is a script for correcting repetitive elements in pacbio assemblies with Illumina data

fix_repeats_ill.py
This script maps (single-end) Illumina reads back to repetitive regions with no Illumina coverage
and corrects the errors found.

USAGE: python fix_repeats_ill.py -c <coverage.txt> -g <genome.fa> -r <reads.fq> -w <working_dir> -o <out_file>
Where coverage.txt is a bed file of the coverage at all bases
(can be generated using genomeCoverageBed -d -ibam aln.sorted.bam -g ref.fasta > coverage.txt)
genome.fa is the reference to be corrected
reads.fq are the illumina reads (can be gzipped)
working_dir is where to put intermediate files
and out_file is the place to write the corrected genome

''', epilog="Thanks for using Contiguity")
parser.add_argument('-c', '--coverage', action='store', help='bed file of coverage at each base')
parser.add_argument('-r', '--read_file', action='store', help='read file (.fastq, .fastq.gz)')
parser.add_argument('-r2', '--read_file_2', default=None, action='store', help='read file (.fastq, .fastq.gz)')
parser.add_argument('-g', '--genome', action='store', help='FASTA file of genome to be corrected')
parser.add_argument('-w', '--working_dir', action='store', help='working directory')
parser.add_argument('-o', '--output_fasta', action='store', help='place to write corrected FASTA')


args = parser.parse_args()


try:
    os.makedirs(args.working_dir)
except:
    pass

correct_regions(args.genome, args.read_file, args.coverage, args.working_dir, args.output_fasta, args.read_file_2)
