#!/usr/bin/env python

import subprocess
import sys
import os

def reorientate_fasta(fasta_file, log_file, out_file, working_dir):
    with open(log_file) as f:
        if f.readline() == 'all_good':
            need_reorient = False
        else:
            need_reorient = True
            subprocess.Popen('makeblastdb -dbtype nucl -in ' + log_file + ' -out ' + working_dir + '/tempdb', shell=True).wait()
            subprocess.Popen('blastn -query ' + fasta_file + ' -db ' + working_dir + '/tempdb -outfmt 6 -out ' + working_dir + '/blast_temp.out', shell=True).wait()
    seqDict = {}
    first = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                seqDict[name] = ''
                first[name] = None
            else:
                seqDict[name] += line.rstrip()
    if need_reorient:
        with open(working_dir + '/blast_temp.out') as f:
            for line in f:
                query, subject, length, ident, mm, indel, qstart, qstop, rstart, rstop, evalue, bitscore = line.split()
                if query + '|quiver' == subject and int(qstart) < int(qstop) and int(rstart) < int(rstop) and int(rstart) < 10 and float(ident) > 90:
                    if first[query] is None:
                        first[query] = map(int, (qstart, qstop, rstart, rstop))
    out = open(out_file, 'w')
    with open(log_file) as f:
        for line in f:
            if line.startswith('>') and not line.rstrip()[1:] in first:
                sys.sterr.write('No hit found for contig ' + line.rstrip()[1:] + '\n')
                sys.exit()
    for i in seqDict:
        if first[i] is None:
            newstart = 0
        else:
            newstart = first[i][0] - first[i][2]
        out.write('>' + i[:8] + 'p' + i[9:-7] + '\n')
        seq = seqDict[i][newstart:] + seqDict[i][:newstart]
        for k in range(0, len(seq), 80):
            out.write(seq[k:k+80] + '\n')
    out.close()

fasta_file = sys.argv[1]
log_file = sys.argv[2]
out_file = sys.argv[3]
working_dir = sys.argv[4]


if not os.path.exists(working_dir):
    os.makedirs(working_dir)

reorientate_fasta(fasta_file, log_file, out_file, working_dir)
