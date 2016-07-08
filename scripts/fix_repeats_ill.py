import argparse
import numpy
import sys
import subprocess
import os


def get_low_coverage(coverage_file):
    with open(coverage_file) as cov:
        cov_dict = {}
        cov_vals = []
        for line in cov:
            ref, pos, cov = line.split()
            if not ref in cov_dict:
                cov_dict[ref] = []
            cov_dict[ref].append(int(cov))
            cov_vals.append(int(cov))
        average = numpy.average(cov_vals)
        std = numpy.std(cov_vals)
        low_cov = {}
        cov_cutoff = max([5, average - 3 * std])
        print cov_cutoff
        for i in cov_dict:
            low_cov[i] = []
            minval = None
            maxval = None
            for pos, j in enumerate(cov_dict[i]):
                if j < cov_cutoff:
                    if minval is None:
                        minval = pos - 100
                    maxval = pos
                else:
                    if not maxval is None and pos > maxval + 100:
                        low_cov[i].append((minval, maxval + 100))
                        minval = None
                        maxval = None
    return low_cov



def correct_regions(fasta_file, read_file, low_cov, working_dir, out_file):
    with open(fasta_file) as fasta:
        seqDict = {}
        for line in fasta:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                seqDict[name] = ''
            else:
                seqDict[name] += line.rstrip()
    out = open(out_file, 'w')
    for i in low_cov:
        for j in low_cov[i]:
            old_seq = seqDict[i][j[0]:j[1]]
            with open(working_dir + '/ref.fa', 'w') as ref:
               ref.write('>ref\n' + old_seq)
            subprocess.Popen('bwa index ' + working_dir + '/ref.fa', shell=True).wait()
            subprocess.Popen('bwa mem -t 4 ' + working_dir + '/ref.fa ' + read_file + ' > ' + working_dir + '/ref.aln.sam', shell=True).wait()
            subprocess.Popen('samtools view -bS ' + working_dir + '/ref.aln.sam > ' + working_dir + '/ref.aln.bam', shell=True).wait()
            subprocess.Popen('samtools sort -o ' + working_dir + '/ref.sort.bam ' + working_dir + '/ref.aln.bam', shell=True).wait()
            subprocess.Popen('samtools mpileup -L100000 -d100000 -uf "' + working_dir + '/ref.fa" ' + working_dir + \
                            '/ref.sort.bam | bcftools call -cv -Ob > "' + working_dir + '/ref.bcf"', shell=True).wait()
            subprocess.Popen('bcftools view "' + working_dir + '/ref.bcf" | vcfutils.pl varFilter > "' + working_dir + '/ref.vcf"', shell=True).wait()
            subprocess.Popen('bgzip -c "' + working_dir + '/ref.vcf" > "' + working_dir + '/ref.vcf.gz"', shell=True).wait()
            subprocess.Popen('tabix -p vcf "' + working_dir + '/ref.vcf.gz"', shell=True).wait()
            subprocess.Popen('cat "' + working_dir + '/ref.fa" | vcf-consensus "' + working_dir + '/ref.vcf.gz" > "' + working_dir + '/new_ref.fasta"', shell=True).wait()
            seq = ''
            with open(working_dir + '/new_ref.fasta') as fasta:
                for line in fasta:
                    if not line.startswith('>'):
                        seq += line.rstrip()
            full_seq = seqDict[i]
            seqDict[i] = full_seq[:j[0]] + seq + full_seq[j[1]:]
        out.write('>' +2 i + '\n')
        seq = seqDict[i]
        for j in range(0, len(seq), 80):
            out.write(seq[j:j+80] + '\n')
    out.close()



try:
    coverage_file, ref_file, read_file, working_dir, out_file = sys.argv[1:]
except:
    print '''
fix_repeats_ill.py
This script maps (single-end) Illumina reads back to repetitive regions with no Illumina coverage
and corrects the errors found.

USAGE: python fix_repeats_ill.py <coverage.txt> <genome.fa> <reads.fq> <working_dir> <out_file>
Where coverage.txt is a bed file of the coverage at all bases
(can be generated using genomeCoverageBed -d -ibam aln.sorted.bam -g ref.fasta > coverage.txt)
genome.fa is the reference to be corrected
reads.fq are the illumina reads (can be gzipped)
working_dir is where to put intermediate files
and out_file is the place to write the corrected genome
'''
low_cov = get_low_coverage(coverage_file)
try:
    os.makedirs(working_dir)
except:
    pass

correct_regions(ref_file, read_file, low_cov, working_dir, out_file)