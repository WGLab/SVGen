# insert_SVs.py
# C: Sep 29, 2015
# M: Dec  8, 2015
# A: Leandro Lima <leandrol@usc.edu>


import os, sys, argparse, linecache
from random import randint
from operator import itemgetter

prog_name = 'insert_SVs.py'


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': ''}



def reverse_complement(seq):
    return ''.join(complement[base] for base in reversed(seq))



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
    return seq[:start] + reverse_complement(seq[start:end+1:]) + seq[end+1:]



def insert_trans(large_seq, small_seq, start, zero_based=False):
    if not zero_based:
        start -= 1
    return large_seq[:start] + small_seq + large_seq[start:]



def calculate_positions(start, end, line_len):
    '''This function calculate the positions of starting line, starting position, ending line and ending position,
    given starting and ending position in a chromosome, assuming that positions are 1-based.
    '''
    start_line = ((start-1) / line_len) + 1
    end_line   = ((end-1)   / line_len) + 1
    start_pos  = start - line_len*(start_line - 1) - 1
    end_pos    = end   - line_len*(end_line   - 1) - 1 
    return start_line, start_pos, end_line, end_pos



def get_subseq_from_fasta(fasta_input_path, start_line, start_pos, end_line, end_pos):
    if start_line < 1:
        print 'Error! Start line has to be greater than 1.'
        return None
    if end_line < start_line:
        print 'Error! End line has to be greater than start line.'
        return None
    line_number = start_line
    seq = linecache.getline(fasta_input_path, line_number)[start_pos:].strip()
    line_number += 1
    while line_number < end_line:
        seq += linecache.getline(fasta_input_path, line_number).strip()
        line_number += 1
    if start_line == end_line:
        seq = linecache.getline(fasta_input_path, start_line)[start_pos:end_pos+1].strip()
        linecache.clearcache()
        return seq
    seq += linecache.getline(fasta_input_path, line_number)[:end_pos+1].strip()
    linecache.clearcache()
    return seq



def write_fasta(chrom_name, output_fasta, seq, size_per_line, fasta_label=None):
    if fasta_label is None:
        output_fasta.write('>' + chrom_name + '\n')
    else:
        output_fasta.write('>' + fasta_label + '\n')
    i = 0
    while i < len(seq):
        output_fasta.write(seq[i:i+size_per_line] + '\n')
        i += size_per_line
    output_fasta.close()



def main():

    parser = argparse.ArgumentParser(description='This program inserts structural variants from a BED file into a FASTA file.', prog=prog_name)
    parser.add_argument('--fasta_input', '-i',  required=True, metavar='input.fasta',  type=file, help='Fasta file to be changed with SVs.')
    parser.add_argument('--fasta_output', '-o', required=True, metavar='output.fasta', type=argparse.FileType('w'), help='Fasta file to be created with SVs.')
    parser.add_argument('--bed', dest='sv_bed', required=True, metavar='SVs.bed', type=file, help='BED file with SVs to be inserted.')
    parser.add_argument('--chrom_lens', required=True, type=file, metavar='chrom_lengths_file', dest='chrom_lens_file', help='Text file with chromosome lengths.')
    # parser.add_argument('--fasta_type', '-f',  required=True, metavar='fasta_type',  type=str, help='A [single] fasta file or [multiple]?')
    parser.add_argument('--chrom', required=True, type=str, metavar='chromosome_name', dest='chromosome_name', help='Chromosome.')
    parser.add_argument('--fasta_label', required=False, type=str, metavar='fasta_label', dest='fasta_label', help='Name to label fasta sequence.')

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()


    chrom_lens = {}
    with args.chrom_lens_file as chrom_lens_file_lines:
        for line in chrom_lens_file_lines:
            try:
                chrom, chrom_len = line.split()
                chrom_lens[chrom] = int(chrom_len)
            except:
                pass


    lines = args.fasta_input.read().split('\n')
    while lines[-1] == '':
        lines.pop()

    if not lines[0].startswith('>'):
        print 'Invalid fasta file!!'
        sys.exit(1)
    else:
        chrom_label = lines[0][1:]


    chrom_name = args.chromosome_name #lines[0].replace('>', '')

    fasta_size_per_line = len(lines[1])

    # if args.fasta_type == 'multiple':
    
    chrom_seq = ''
    LEN_CHROM = 0
    for line in lines[1:]:
        LEN_CHROM += len(line)
        chrom_seq += line.upper()
    # chrom_seq += line

    lines = args.sv_bed.read().split('\n')
    while lines[-1] == '':
        lines.pop()

    sv_regions = []
    for line in lines:
        elems = line.split()
        chrom, start, end, sv_type = elems[:4]
        if len(elems) > 4:
            translocation_info = elems[4]
            sv_regions.append([chrom, int(start), int(end), sv_type, translocation_info])
        else:
            sv_regions.append([chrom, int(start), int(end), sv_type])

    sv_regions = sorted(sv_regions, key=itemgetter(1), reverse=True)
    # print sv_regions

    print 'Inserting SVs in sequence.'
    for sv_region in sv_regions:
        chrom, start, end, sv_type = sv_region[:4]
        if chrom_name == chrom:
            if args.verbose:
                print 'Inserting %s in %s:%s-%s.' % (sv_type.upper(), chrom, start, end)
            if len(sv_region) > 4:
                translocation_info = sv_region[4]
            if sv_type.upper().startswith('DEL'):
                chrom_seq = insert_del(chrom_seq, start, end)
            elif sv_type.upper().startswith('DUP'):
                chrom_seq = insert_dup(chrom_seq, start, end)
            elif sv_type.upper().startswith('INV'):
                chrom_seq = insert_inv(chrom_seq, start, end)
            elif sv_type.upper().startswith('UNBTR') or sv_type.upper().startswith('BALTR'):
                chrom_tr = translocation_info.split(':')[0]
                start_tr, end_tr = translocation_info.split(':')[1].split('-')
                start_line, start_pos, end_line, end_pos = calculate_positions(start, end, fasta_size_per_line)
                trans_chrom_filename  = args.fasta_input.name.replace('chr'+chrom_name, 'chr'+chrom_tr)
                # trans_chrom_filename2 = args.fasta_input.name.replace(chrom_name, 'chr'+chrom_tr)
                if not os.path.isfile(trans_chrom_filename): # and not os.path.isfile(trans_chrom_filename2):
                    # print 'Error! I was not able to find files [%s|%s].' % (trans_chrom_filename, trans_chrom_filename2)
                    print 'Error! I was not able to find files [%s].' % (trans_chrom_filename)
                    sys.exit(0)
                # elif not os.path.isfile(trans_chrom_filename):
                    # trans_chrom_filename = trans_chrom_filename2
                if args.verbose:
                    print 'Getting sequence %s from file [%s].' % (translocation_info, trans_chrom_filename)
                trans_seq = get_subseq_from_fasta(trans_chrom_filename, start_line, start_pos, end_line, end_pos)
                if sv_type.upper().startswith('BALTR'):
                    chrom_seq = insert_del(chrom_seq, start, end)
                chrom_seq = insert_trans(chrom_seq, trans_seq, start)


    print 'Writing to output file [%s].' % args.fasta_output.name
    write_fasta(chrom_name, args.fasta_output, chrom_seq, fasta_size_per_line, chrom_label)
    print 'Simulations for', chrom_name, 'done.'
    print 'Total size:', LEN_CHROM
    print 


if __name__ == '__main__':
    main()


