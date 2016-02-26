#!/usr/bin/env python
# checkCirclator.py   Written by: Mitchell Sullivan   mjsull@gmail.com
# Supervisor: Dr. Harm Van Bakel
# Version 0.0.1 2016.01.12
# License: GPLv3

import sys



def checkLog(circ_direct, outname, seqlog, assembly_no):
    min_dist = 5000
    seqDict = {}
    name = None
    with open(circ_direct + '/06.fixstart.fasta') as f: # get the final output of circlator
        for line in f:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                seqDict[name] = ''
            else:
                seqDict[name] += line.rstrip()
    seqDictOri = {}
    with open(circ_direct + '/00.input_assembly.fasta') as f: # get the initial input for circlator - this is done so we can readd unitigs thrown out by circlator
        for line in f:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                seqDictOri[name] = ''
            else:
                seqDictOri[name] += line.rstrip()
    circ_set = set()
    with open(circ_direct + '/04.merge.circularise.log') as f: # determine whether circlator has circularised the unitig
        first = True
        for line in f:
            if first:
                first = False
            else:
                mc1, mc2, name, one, two, three, circed = line.split()
                if circed == '1':
                    circ_set.add(name)
    reorient_dict = set()
    with open(circ_direct + '/06.fixstart.log') as f: # determine whether circlator has reorientated the unitig
        first = True
        for line in f:
            if first:
                first = False
            else:
                cbf1, cbf2, cbf3, name, break_point, gn, gr, nn, skipped = line.split()
                if break_point == '-':
                    pass
                else:
                    reorient_dict[name] = int(break_point)
    out = open(outname, 'w')
    out_log = open(seqlog, 'w')
    at_least_one = False
    for i in seqDictOri: # build new name for contig
        contig_name = '>u' + i.split('|')[0].split('_')[1].zfill(5)
        if i in circ_set:
            contig_name += 'c'
        else:
            contig_name += 'x'
        if i in reorient_dict:
            contig_name += 'r'
        else:
            contig_name += 'x'
        contig_name += 'xx_'
        if not i in seqDict:
            contig_name += 'g'
        elif len(seqDictOri[i]) >= 1000000 and i in circ_set:
            contig_name += 'c'
        elif i in circ_set:
            contig_name += 'p'
        else:
            contig_name += 'o'
        if not i in seqDict:
            seq = seqDictOri[i]
        elif not i in circ_set:
            seq = seqDict[i]
        elif not i in reorient_dict:
            seq = seqDict[i][len(seq)/2:] + seqDict[i][:len(seq)/2]
        elif i in reorient_dict and reorient_dict[i] >= min_dist and reorient_dict[i] <= len(seqDict[i]) - min_dist:
            seq = seqDict[i]
        elif i in reorient_dict and reorient_dict[i] >= len(seqDict[i]) /4 and reorient_dict[i] <= len(seqDict[i]) * 3 / 4:
            seq = seqDict[i]
        else:
            at_least_one = True
            out_log.write(contig_name + '\n')
            out_log.write(seqDict[i])
            seq = seqDict[i][len(seq)/2:] + seqDict[i][:len(seq)/2]
        contig_name += '_' + assembly_no
        out.write(contig_name + '\n')
        for j in range(0, len(seq), 60):
            out.write(seq[j:j+60] + '\n')
    if not at_least_one: # if no contigs needed to be reorientated write all good to file - post quiver will check this and then not bother running BLAST to change files back
        out_log.write('all_good')
    out_log.close()
    out.close()

if len(sys.argv) == 5:
    checkLog(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
else:
    print '''checkCirclator.py
checkCirclator.py checks circlator output to see if reorientation point was chosen near a contig start or end
if it was circlator reorientates contig so that quiver may correct poor quality ends

After Quiver has been run - postquiver_orient_correct.py quiver_out.fasta outfile2.fasta
may be run to reorientate contigs and correct names.

USAGE: checkCirclator.py circlator_output_dir outfile.fasta outfile2.fasta'''
