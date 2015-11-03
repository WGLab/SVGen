# insert_SVs.py
# C: Sep 29, 2015
# M: Nov  2, 2015
# A: Leandro Lima <leandrol@usc.edu>


import sys
from random import randint
from operator import itemgetter


def insert_del(seq, start, end, zero_based=False):
    if end <= start:
        print 'Error! End position must be greater than Start.'
        return seq
    if not zero_based:
        start -= 1
        end   -= 1
    return seq[:start] + seq[end+1:]



def insert_dup(seq, start, end, zero_based=False):
    if end <= start:
        print 'Error! End position must be greater than Start.'
        return seq
    if not zero_based:
        start -= 1
        end   -= 1
    return seq[:end+1] + seq[start:]



def insert_inv(seq, start, end, zero_based=False):
    if end <= start:
        print 'Error! End position must be greater than Start.'
        return seq
    if not zero_based:
        start -= 1
        end   -= 1
    return seq[:start] + seq[start:end+1:][::-1] + seq[end+1:]
 


def write_fasta(chrom_name, output_fasta_name, seq, size_per_line, variant_type):
    chrom_fasta_file = open(output_fasta_name, 'w')
    chrom_fasta_file.write('>' + chrom_name + '_' + variant_type + '\n')
    i = 0
    while i < len(seq):
        chrom_fasta_file.write(seq[i:i+size_per_line] + '\n')
        i += size_per_line
    chrom_fasta_file.close()



def main():

    if not len(sys.argv) == 4:
        print 'Usage: python simulate_CNVs.py [fasta_input] [fasta_output] [bed_file]'

    fasta_input      = sys.argv[1]
    fasta_output     = sys.argv[2]
    sv_bed_filename  = sys.argv[3]

    lines = open(fasta_input).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    if not lines[0].startswith('>'):
        print 'Invalid fasta file!!'
        sys.exit(1)

    chrom_name = lines[0].replace('>', '')

    fasta_size_per_line = len(lines[1])

    chrom_seq = ''
    LEN_CHROM = 0
    for line in lines[1:]:
        LEN_CHROM += len(line)
        chrom_seq += line.upper()

    lines = open(sv_bed_filename).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    cnv_regions = []
    for line in lines:
        chrom, start, end, cnv_type = line.split()
        cnv_regions.append([chrom, int(start), int(end), cnv_type])

    cnv_regions = sorted(cnv_regions, key=itemgetter(1), reverse=True)

    print 'Inserting CNVs in sequence.'
    for cnv_region in cnv_regions:
        chrom, start, end, cnv_type = cnv_region
        if cnv_type.upper().startswith('DEL'):
            chrom_seq = insert_del(chrom_seq, start, end)
        elif cnv_type.upper().startswith('DUP'):
            chrom_seq = insert_dup(chrom_seq, start, end)
        elif cnv_type.upper().startswith('INV'):
            chrom_seq = insert_inv(chrom_seq, start, end)


    print 'Writing to output file [%s].' % fasta_output
    write_fasta(chrom_name, fasta_output, chrom_seq, fasta_size_per_line, 'CNVs')
    print 'Simulations for', chrom_name, 'done.'
    print 'Total size:', LEN_CHROM
    print 


if __name__ == '__main__':
    main()


