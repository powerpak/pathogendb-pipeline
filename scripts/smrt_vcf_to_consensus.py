import sys

try:
    min_qual = float(sys.argv[3])
except:
    sys.sterr.write('USAGE smrt_vcf_to_consensus.py input.fasta input.vcf min_quality\n')
    sys.exit()

seq_dict = {}
change_dict = {}
order_list = []
with open(sys.argv[1]) as fasta:
    for line in fasta:
        if line.startswith('>'):
            name = line.rstrip()[1:]
            seq_dict[name] = ''
            change_dict[name] = []
            order_list.append(name)
        else:
            seq_dict[name] += line.rstrip()


with open(sys.argv[2]) as vcf:
    for line in vcf:
        if not line.startswith('#'):
            chrom, pos, id, ref, alt, qual, filt, inf = line.split()
            if float(qual) > min_qual:
                change_dict[chrom].append((pos, ref, alt))

for i in order_list:
    var_list = change_dict[i]
    var_list.sort()
    seq = seq_dict[i]
    mod = 0
    for j in var_list:
        pos, ref, alt = j
        pos = int(pos) + mod
        if alt[0] == 'I':
            seq = seq[:pos] + alt[1:] + seq[pos:]
            mod += len(alt) - 1
        elif alt[0] == 'D':
            bases = int(alt[1:])
            seq = seq[:pos] + seq[pos+bases:]
            mod -= bases
        else:
            seq = seq[pos-1:] + alt + seq[pos-1 + len(alt):]
    sys.stdout.write('>' + i + '\n')
    for j in range(0, len(seq), 80):
        sys.stdout.write(seq[j:j+80] + '\n')


