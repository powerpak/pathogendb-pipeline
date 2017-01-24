#!/usr/bin/env python

import sys
import argparse
import subprocess

def reverse_compliment(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def trim_contigs(args):
    seqDict = {}
    with open(args.fasta) as fasta:
        for line in fasta:
            if line.startswith('>'):
                name = line.split()[0][1:]
                seqDict[name] = ''
            else:
                seqDict[name] += line.rstrip()
    if not args.trim is None:
        for num in range(0, len(args.trim), 2):
            trim_name = args.trim[num]
            trim = int(args.trim[num+1])
            for i in seqDict:
                if trim_name in i:
                    sys.stdout.write('Trimmed contig: ' + i + '\n')
                    if trim > 0:
                        seqDict[i] = seqDict[i][trim:]
                    else:
                        seqDict[i] = seqDict[i][:trim]
    with open('break.fasta', 'w') as f:
        for i in seqDict:
            f.write('>' + i + '\n' + seqDict[i] + '\n')
    subprocess.Popen('makeblastdb -dbtype nucl -in break.fasta -out breakdb', shell=True, stdout=subprocess.PIPE).wait()
    subprocess.Popen('blastn -query break.fasta -db breakdb -outfmt 6 -out break.out', shell=True).wait()
    min_length = args.min_length
    min_ident = args.min_ident
    overlap_dict = {}
    wobble = args.wobble
    merge_list = []
    with open('break.out') as blast:
        for line in blast:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length = map(int, [qstart, qstop, rstart, rstop, length])
            eval, bitscore, ident = map(float, [eval, bitscore, ident])
            if qstart <= wobble and length >= min_length and ident >= min_ident and qstart != rstart and rstop >= len(seqDict[query]) - wobble and query == subject:
                print line.rstrip()
                if query in overlap_dict:
                    if qstop > overlap_dict[query]:
                        overlap_dict[query] = qstop
                else:
                    overlap_dict[query] = qstop
            elif qstart <= wobble and length >= min_length and ident >= min_ident and rstop >= len(seqDict[subject]) - wobble and query != subject:
                merge_list.append([subject, '+', query, '+', qstop])
            elif rstart <= wobble and length >= min_length and ident >= min_ident and qstop >= len(seqDict[query]) - wobble and query != subject:
                merge_list.append([query, '+', subject, '+', rstop])
            elif qstart <= wobble and length >= min_length and ident >= min_ident and rstop <= wobble and query != subject:
                merge_list.append([subject, '-', query, '+', qstop])
            elif qstop >= len(seqDict[query]) - wobble and length >= min_length and ident >= min_ident and rstart >= len(seqDict[subject]) - wobble:
                merge_list.append([query, '+', subject, '-', rstart - rstop])
    shorten_names = {}
    for i in seqDict:
        if i.startswith('unitig_'):
            shorten_names[i] = i.split('|')[0].split('_')[1]
        elif i.startswith('u'):
            shorten_names[i] = str(int(i[:6]))
    if not args.merge is None:
        for i in range(0, len(args.merge), 4):
            count = i / 4 + 1
            contig1, dir1, contig2, dir2 = args.merge[i:i+4]
            for j in merge_list:
                merge1, mdir1, merge2, mdir2, trimit = j
                if (contig1 == merge1 or contig1 == shorten_names[merge1]) and dir1 == mdir1 and (contig2 == merge2 or contig2 == shorten_names[merge2]) and dir2 == mdir2:
                    mergeit = j
                elif (contig2 == merge1 or contig2 == shorten_names[merge1]) and dir1 != mdir1 and (contig1 == merge2 or contig1 == shorten_names[merge2]) and dir2 == mdir2:
                    mergeit = j
            merge1, mdir1, merge2, mdir2, trimit = mergeit
            new_contig_name = 'merged_' + str(count)
            seq1 = seqDict[merge1]
            if mdir1 == '-':
                seq1 = reverse_compliment(seq1)
            seq2 = seqDict[merge2]
            if mdir2 == '-':
                seq2 = reverse_compliment(seq2)
            new_contig_seq = seq1 + seq2[trimit:]
            del seqDict[merge1]
            del seqDict[merge2]
            sys.stdout.write('merged ' + str(merge1) + ' ' + str(merge2) + ' overlap:' + str(trimit) + '\n')
            seqDict[new_contig_name] = new_contig_seq
            shorten_names[new_contig_name] = new_contig_name
    with open(args.output, 'w') as outfile:
        for i in seqDict:
            if args.remove is None or (not i in args.remove and not shorten_names[i] in args.remove):
                if i in overlap_dict:
                    seq = seqDict[i][overlap_dict[i]:]
                    sys.stdout.write(i + ' length:' + str(len(seq)) + ' trimmed:' + str(overlap_dict[i]) + '\n')
                else:
                    seq = seqDict[i]
                    sys.stdout.write(i + ' length:' + str(len(seq)) + ' trimmed:0\n')
                outfile.write('>' + i + '\n')
                for j in range(0, len(seq), 60):
                    outfile.write(seq[j:j+60] + '\n')
            else:
                sys.stdout.write('filtered: ' + i + '\n')


parser = argparse.ArgumentParser(description='''
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
parser.add_argument("-o", "--output", help="prefix used for output files", required=True)
parser.add_argument("-f", "--fasta", help="FASTA of genome.", metavar="genome.fasta", required=True)
parser.add_argument("-w", "--wobble", type=int, default=50, help="Minimum length of alignment to count as repeat.")
parser.add_argument("-i", "--min_ident", type=float, default=96.0, help="Routes PHASTER through HTTP proxy.")
parser.add_argument("-l", "--min_length", type=int, default=500, help="path to database of phage proteins")
parser.add_argument("-r", "--remove", nargs="+", help="unitigs to filter")
parser.add_argument("-m", "--merge", nargs="+", help="merge these unitigs - e.g. unitig1 + unitig3 - will merge the positive strand of unitig 1 and negative strand of unitig3 if overlap is found")
parser.add_argument("-t", "--trim", nargs="+", help="Trim edge of contigs before finding overlap")
args = parser.parse_args()

trim_contigs(args)




