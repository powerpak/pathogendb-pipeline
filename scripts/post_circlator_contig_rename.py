#!/usr/bin/env python
# checkCirclator.py   Written by: Mitchell Sullivan   mjsull@gmail.com
# Supervisor: Dr. Harm Van Bakel
# Version 0.0.1 2016.01.12
# License: GPLv3

import sys
import os


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
    if os.path.exists(circ_direct + '/00.input_assembly.fasta'):
        with open(circ_direct + '/00.input_assembly.fasta') as f: # get the initial input for circlator - this is done so we can readd unitigs thrown out by circlator
            for line in f:
                if line.startswith('>'):
                    name = line.rstrip()[1:]
                    seqDictOri[name] = ''
                else:
                    seqDictOri[name] += line.rstrip()
        curated = False
    else:
        seqDictOri = seqDict
        curated = True
        name_list = list(seqDict)
        name_list.sort(key=lambda x: len(seqDict[x]), reverse=True)
        name_dict = {}
        for num, i in enumerate(name_list):
            name_dict[i] = '>u' + str(num).zfill(5)
    circ_set = set()
    if curated:
        circ_set = set(seqDict)
    else:
        with open(circ_direct + '/04.merge.circularise.log') as f: # determine whether circlator has circularised the unitig
            first = True
            for line in f:
                if first:
                    first = False
                else:
                    mc1, mc2, name, one, two, three, circed = line.split()
                    if circed == '1':
                        circ_set.add(name)
    reorient_dict = {}
    with open(circ_direct + '/06.fixstart.log') as f: # determine whether circlator has reorientated the unitig, keep track of breakpoint
        first = True
        for line in f:
            if first:
                first = False
            else:
                name, break_point, gn, gr, nn, skipped = line.split(']')[1].split()
                if skipped ==   'skipped':
                    pass
                else:
                    if break_point == '-':
                        break_point = 0
                    reorient_dict[name] = int(break_point)
    out = open(outname, 'w')
    out_log = open(seqlog, 'w')
    at_least_one = False
    out_list = []
    for q in seqDictOri: # build new name for contig
        gotit = False
        merged = False
        for j in seqDict:
            if j.startswith(q):
                i = j
                gotit = True
                break
            elif q in j:
                merged = True
        if not gotit:
            i = q
        if curated:
            contig_name = name_dict[i]
            manual = True
        elif i.startswith('unitig'):
            contig_name = '>u' + i.split('|')[0].split('_')[1].zfill(5)
            manual = False
        else:
            contig_name = '>' + i[:6]
            manual = True
        if i in circ_set:
            contig_name += 'c'
        else:
            contig_name += 'x'
        if i in reorient_dict:
            contig_name += 'r'
        else:
            contig_name += 'x'
        if manual:
            contig_name += 'xm_'
        else:
            contig_name += 'xx_'
        if not i in seqDict:
            if merged:
                contig_name += 'm' # label contig as merged into larger contig
            else:
                contig_name += 'g' # recover tossed out contigs (and label them as "g" for garbage)
        elif i in seqDict and len(seqDict[i]) >= 1000000 and i in circ_set: # assign contig to class (circularised and large = chromosome, circularised and small = plasmid, not circularised = other)
            contig_name += 'c'
        elif i in circ_set:
            contig_name += 'p'
        else:
            contig_name += 'o'
        if not i in seqDict:
            seq = seqDictOri[i]
        elif not i in circ_set: # if contig was not circularised do not reorientate
            seq = seqDict[i]
        elif not i in reorient_dict: # if contig was not reorientated, re-reorientate for quiver
            seq = seqDict[i][len(seq)/2:] + seqDict[i][:len(seq)/2]
        # if the contig was reorientated far enough away from the initial contig start do not re-reorientate
        elif i in reorient_dict and reorient_dict[i] >= min_dist and reorient_dict[i] <= len(seqDict[i]) - min_dist:
            seq = seqDict[i]
        # if the contig is short and was reorientated roughly in the middle do no re-reorientate
        elif i in reorient_dict and reorient_dict[i] >= len(seqDict[i]) /4 and reorient_dict[i] <= len(seqDict[i]) * 3 / 4:
            seq = seqDict[i]
        else: # If the contig has been reorientated near the start, re-reorienate and keep a record of inital position so that it can be orientated back to DNAA after quiver
            at_least_one = True
            out_log.write(contig_name + '\n')
            out_log.write(seqDict[i])
            seq = seqDict[i][len(seq)/2:] + seqDict[i][:len(seq)/2]
        contig_name += '_' + assembly_no
        out_list.append((contig_name, seq))
    out_list.sort()
    for contig_name, seq in out_list:
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
