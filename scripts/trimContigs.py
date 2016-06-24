#!/usr/bin/env python

import sys
import subprocess



def trim_contigs(in_fasta, out_fasta):
    seqDict = {}
    with open(in_fasta) as fasta:
        for line in fasta:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                seqDict[name] = ''
            else:
                seqDict[name] += line.rstrip()
    subprocess.Popen('makeblastdb -dbtype nucl -in ' + sys.argv[1] + ' -out breakdb', shell=True, stdout=subprocess.PIPE).wait()
    subprocess.Popen('blastn -query ' + sys.argv[1] + ' -db breakdb -outfmt 6 -out break.out', shell=True).wait()
    min_length = 500
    min_ident = 95
    overlap_dict = {}
    wobble = 5
    with open('break.out') as blast:
        for line in blast:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length = map(int, [qstart, qstop, rstart, rstop, length])
            eval, bitscore, ident = map(float, [eval, bitscore, ident])
            if qstart <= wobble and length >= min_length and ident >= min_ident and qstart != rstart and rstop >= len(seqDict[query]) - wobble and query == subject:
                if query in overlap_dict:
                    if qstart > overlap_dict[query]:
                        overlap_dict[query] = qstop
                else:
                    overlap_dict[query] = qstop
    filterset = set()
    for i in sys.argv[3:]:
        filterset.add('u' + i.zfill(5))
    with open(out_fasta, 'w') as outfile:
        for i in seqDict:
            if not i[:6] in filterset:
                if i in overlap_dict:
                    print i, overlap_dict[i]
                    seq = seqDict[i][overlap_dict[i]:]
                else:
                    print i, 0
                    seq = seqDict[i]
                outfile.write('>' + i + '\n')
                for j in range(0, len(seq), 60):
                    outfile.write(seq[j:j+60] + '\n')


if len(sys.argv) < 3:
    sys.stdout.write('''
trimContigs.py
writtne by Mitchell Sullivan
mjsull@gmail.com for help

USAGE: trimContigs.py input_fasta.fa output_fasta.fa <contig_numbers>

EXAMPLE A: trimContigs.py input_fasta.fa output_fasta.fa
If the edges of input_fasta.fa overlap trim overlap from start of the contig or else do nothing

EXAMPLE B: trimContigs.py input_fasta.fa output_fasta.fa 3 4
If the edges of input_fasta.fa overlap trim overlap from start of the contig
Additionally remove u00003 and u00004 from the assembly

REQUIRES BLAST BE INSTALLED ON THE COMMAND LINE

WILL OVERWRITE breakdb.out breakdb.nhr breakdb.nin breakdb.nsq - do not run concurrently from the same directory.

''')
else:
    trim_contigs(sys.argv[1], sys.argv[2])




