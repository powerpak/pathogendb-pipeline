#!/usr/bin/env python

import requests
import json
import time
import subprocess
import random
import argparse
import os

# takes HSL value and converts to RGB format for bed files
def hsl_to_colorstr(h, s, l):
    c = (1 - abs(2*l - 1)) * s
    x = c * (1 - abs(h *1.0 / 60 % 2 - 1))
    m = l - c/2
    if h < 60:
        r, g, b = c + m, x + m, 0 + m
    elif h < 120:
        r, g, b = x + m, c+ m, 0 + m
    elif h < 180:
        r, g, b = 0 + m, c + m, x + m
    elif h < 240:
        r, g, b, = 0 + m, x + m, c + m
    elif h < 300:
        r, g, b, = x + m, 0 + m, c + m
    else:
        r, g, b, = c + m, 0 + m, x + m
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    return str(r) + ',' + str(g) + ',' + str(b)


# submit a job to PHASTER and then return the job id
def submit_phaster(fasta, proxy):
    job_ids = []
    with open(fasta, 'r') as f:
        entry = f.readline()
        outlist = []
        if entry.startswith('>'):
            name = entry.rstrip()[1:]
            for line in f:
                if line.startswith('>'):
                    outlist.append((name, entry))
                    name = line.rstrip()[1:]
                    entry = line
                else:
                    entry += line
        elif entry.startswith('LOCUS'):
            name = entry.split()[1]
            for line in f:
                if line.startswith('LOCUS'):
                    outlist.append((name, entry))
                    name = line.split()[1]
                    entry = line
                else:
                    entry += line
        outlist.append((name, entry))
    names = []
    for i in outlist:
        print i[0]
        print i[1][:30]
        if proxy is None:
            r = requests.post('http://phaster.ca/phaster_api', data=i[1])
        else:
            r = requests.post('http://phaster.ca/phaster_api', data=i[1], proxies=proxy)
        jason = json.loads(r.text)
        print jason
        job_ids.append(jason['job_id'])
        names.append(i[0])
    return job_ids, names


# takes a job id and waits for PHASTER job to finish and then writes to bed file
def get_phaster(id, entry_name, outfile, proxy):
    wait_time = 60
    status = 'incomplete'
    while not status == 'Complete':
        time.sleep(wait_time)
        if proxy is None:
            r = requests.get('http://phaster.ca/phaster_api?acc=' + id)
        else:
            r = requests.get('http://phaster.ca/phaster_api?acc=' + id, proxy)
        jason = json.loads(r.text)
        print jason
        status = jason['status']
    out = open(outfile + '.phage.bed', 'a')
    get_phage = False
    for i in jason['summary'].split('\n'):
        print i
        if i.startswith('                                 ---------------------------------'):
            get_phage = True
        elif get_phage:
            try:
                region, region_length, score, specific, phage_start_stop, trna, pn, phpn, hpn, phpp, bpn, att, psn, name, fmcpn, fmcpp, gcp = i.split()
                name = name.split('(')[0]
                start, stop = phage_start_stop.split('-')
                out.write('\t'.join([entry_name, start, stop, name, score, '+']) + '\n')
            except ValueError:
                pass

# Finds repeats in a FASTA file using nucmer
def get_repeats(fasta, prefix, min_ident=85.0, min_length=1000):
    subprocess.Popen('nucmer --maxmatch --nosimplify --prefix=' + prefix + ' ' + fasta + ' ' + fasta, shell=True, stdout=log, stderr=log).wait()
    subprocess.Popen('show-coords -r ' + prefix + '.delta > ' + prefix + '.coords', shell=True, stdout=log, stderr=log).wait()
    repeats = []
    wobble = 5
    with open(prefix + '.coords') as f:
        get_coords = False
        for line in f:
            if line.startswith('=============='):
                get_coords = True
            elif get_coords:
                s1, e1, bar, s2, e2, bar, l1, l2, bar, ident, bar, tag1, tag2 = line.split()
                s1, s2, e1, e2, l1 = map(int, (s1, s2, e1, e2, l1))
                if s2 < e2:
                    orient = '+'
                else:
                    orient = '-'
                    s2, e2, = e2, s2
                ident = float(ident)
                if ident > min_ident and l1 > min_length and not s1 == s2:
                    in_repeat = False
                    for repeat in repeats:
                        s1in, s2in = False, False
                        for i in repeat:
                            if s1 - wobble <= i[0] <= s1 + wobble and e1 - wobble <= i[1] <= e1 + wobble and tag1 == i[2]:
                                s1in = True
                                morient = i[3]
                            if s2 - wobble <= i[0] <= s2 + wobble and e2 - wobble <= i[1] <= e2 + wobble and tag2 == i[2]:
                                s2in = True
                                morient = i[3]
                        if s1in and not s2in:
                            if morient == orient:
                                norient = '+'
                            else:
                                norient = '-'
                            repeat.append((s2, e2, tag2, norient))
                        elif s2in and not s1in:
                            if morient == orient:
                                norient = '+'
                            else:
                                norient = '-'
                            repeat.append((s1, e1, tag1, norient))
                        if s1in or s2in:
                            in_repeat = True
                    if not in_repeat:
                        repeats.append([(s1, e1, tag1, '+'), (s2,e2,tag2, orient)])
    with open(prefix + '.repeats.bed', 'w') as out:
        out.write('track name="Repeats" description="repeat sequence determined by nucmer" itemRgb="On"\n')
        h = 1
        step = 360 / len(repeats)
        for num, i in enumerate(repeats):
            s = 0.4 + 0.2 * random.random()
            l = 0.5 + 0.2 * random.random()
            h += step
            label = 'r' + str(num + 1)
            color = hsl_to_colorstr(h,s,l)
            for j in i:
                out.write('\t'.join([j[2], str(j[0]), str(j[1]), label, '0', j[3], str(j[0]), str(j[1]), color]) + '\n')
        out.close()


# finds islands in a FASTA files using alien_hunter
def get_islands(fasta, prefix):
    subprocess.Popen(os.environ['ALIEN_DIR'] + '/alien_hunter ' + fasta + ' ' + prefix + '.pai', shell=True, stdout=log, stderr=log).wait()
    len_list = []
    name_list = []
    with open(fasta) as fa:
        for line in fa:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                name_list.append(name)
                len_list.append(0)
            else:
                len_list[-1] += len(line.rstrip())
    pai_num = 1
    out = open(prefix + '.pai.bed', 'w')
    out.write('track name="GIs" description="Genomic islands" itemRgb="On"\n')
    with open(prefix + '.pai') as embl:
        for line in embl:
            if line.startswith('FT   misc_feature    '):
                pai_label = 'pai_' + str(pai_num)
                pai_num += 1
                start, stop = line.split()[2].split('..')
                start, stop = int(start), int(stop)
                print start, len_list
                start_pos = 0
                for num, i in enumerate(len_list):
                    if start_pos <= start <= start_pos + i:
                        name = name_list[num]
                        start, stop = str(start - start_pos), str(stop - start_pos)
                        break
                    start_pos += i
                print start,stop
            elif line.startswith('FT                   /colour='):
                color = line.rstrip().split('=')[1].replace(' ', ',')
            elif line.startswith('FT                   /score='):
                score = line.rstrip().split('=')[1]
                out.write('\t'.join([name, start, stop, pai_label, score, '+', start, stop, color]) + '\n')

def get_phage_homebrew(fasta, output, db_path):
    out_list = []
    subprocess.Popen('blastx -outfmt 6 -num_threads 8 -query ' + fasta + ' -db ' + db_path + ' -out ' + output + '.phage.out', shell=True).wait()
    with open(output + '.phage.out') as blast:
        phage_dict = {}
        for line in blast:
            query, subject, ident, length, mm, indel, qstart, qend, rstart, rend, eval, bitscore = line.split()
            if float(ident) >= 33 and int(length) >= 50:
                if not query in phage_dict:
                    phage_dict[query] = set()
                    print qstart, qend, query
                for j in range(min([int(qstart), int(qend)]), max([int(qstart), int(qend)]) + 1):
                    phage_dict[query].add(j)
        gap = 0
        gapmax = 2000
        min_length = 1000
        for i in phage_dict:
            getit = False
            for j in range(0, max(phage_dict[i]) + 1):
                if j in phage_dict[i]:
                    gap = 0
                    if not getit:
                        getit = True
                        phage_start = j
                elif getit:
                    if gap == 0:
                        phage_end = j
                    gap += 1
                    if gap > gapmax:
                        gap = 0
                        getit = False
                        print phage_end, phage_start
                        if phage_end - phage_start >= min_length:
                            out_list.append((i, int(phage_start), int(phage_end)))
            phage_end = j
            if getit:
                if phage_end - phage_start >= min_length:
                    out_list.append((i, int(phage_start), int(phage_end)))
    out = open(output + '.phage.bed', 'w')
    out.write('track name="Phage" description="Phage sequence determined by PHASTER"\n')
    for i in out_list:
        out.write('\t'.join([i[0], str(i[1]), str(i[2]), 'phage', '1', '+']) + '\n')

parser = argparse.ArgumentParser(description='Given a FASTA file will find any or all of phage regions using PHASTER (output to prefix.phage.bed)'
                                             ' repeats found using Nucmer (output to prefix.repeats.bed) and '
                                             'pathogenicity islands using alien_hunter (output to prefix.pai.bed).')
parser.add_argument("-o", "--output_prefix", help="prefix used for output files", required=True)
parser.add_argument("-f", "--fasta", help="FASTA of genome.", metavar="genome.fasta", required=True)
parser.add_argument("-g", "--genbank", help="genbank of genome.", metavar="genome.gbk")
parser.add_argument("-ph", "--phaster", help="Find phage using PHASTER", action='store_true')
parser.add_argument("-p", "--phage", help="Find phage using hombrew phage finder", action='store_true')
parser.add_argument("-i", "--islands",  help="Find genomic islands using alien_hunt", action='store_true')
parser.add_argument("-r", "--repeats", help="Find repeats using nucmer", action='store_true')
parser.add_argument("-m", "--min_ident", default=85.0, type=float, help="Minimum identity of alignment to coun as repeats")
parser.add_argument("-l", "--min_length", default=1000, type=int, help="Minimum length of alignment to count as repeat.")
parser.add_argument("-y", "--proxy", default='False', help="Routes PHASTER through HTTP proxy.")
parser.add_argument("-d", "--phage_db", help="path to database of phage proteins")
args = parser.parse_args()


log = open(args.output_prefix + '.log', 'w')
if args.proxy == 'False':
    proxies = {'http':args.proxy}
else:
    proxies = None
if args.phaster:
    if args.genbank is None:
        job_ids, names = submit_phaster(args.fasta, proxies)
    else:
        job_ids, names = submit_phaster(args.genbank, proxies)
if args.repeats:
    get_repeats(args.fasta, args.output_prefix)
if args.islands:
    get_islands(args.fasta, args.output_prefix)
if args.phage:
    get_phage_homebrew(args.fasta, args.output_prefix, args.phage_db)
if args.phaster:
    out = open(args.output_prefix + '.phage.bed', 'w')
    out.write('track name="Phage" description="Phage sequence determined by PHASTER"\n')
    out.close()
    for i, j in zip(job_ids, names):
        get_phaster(i, j, args.output_prefix, proxies)