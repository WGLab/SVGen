# insert_SVs_and_indels.py
# C: Sep 29, 2015
# M: Apr 26, 2016
# A: Leandro Lima <leandrol@usc.edu>


import os, sys, argparse, linecache
from random import randint
from operator import itemgetter
from copy import copy

prog_name = 'insert_SVs_and_indels.py'


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': ''}

def intervals_overlap(interval1, interval2):
    '''
    Returns a state of overlapping intervals, with the following code:

    0    : no overlap
    1221 : interval 2 is inside interval 1
    2112 : interval 1 is inside interval 2
    1    : any other type of overlap
    '''
    start1, end1 = interval1
    if start1 > end1:
        start1, end1 = end1, start1
    start2, end2 = interval2
    if start2 > end2:
        start2, end2 = end2, start2
    if start1 < start2 and end2 < end1:
        return 1221
    if start2 < start1 and end1 < end2:
        return 2112
    if start1 <= start2 <= end1 or start2 <= start1 <= end2:
        return 1
    return 0



def merge_regions(chrom_regions):
    if not chrom_regions:
        return []
    chrom_regions = sorted(chrom_regions, key = lambda x: x[1])
    result = []
    interval_a = chrom_regions[0] #(a, b)
    new_interval = copy(interval_a) #(a, b)
    interval_set = [interval_a]
    for interval_b in chrom_regions[1:]: # (x, y)
        if interval_b[1] <= new_interval[2]: #x <= b:
            new_interval[2] = max(new_interval[2], interval_b[2]) # b = max(b, y)
            interval_set.append(interval_b)
        else:
            result.append(interval_set) # result.append((a, b))
            interval_a = interval_b # (a, b) = (x, y)
            new_interval = copy(interval_a)
            interval_set = [interval_a]
    result.append(interval_set) # result.append((a, b))
    return result



def reverse_complement(seq):
    return ''.join(complement[base] for base in reversed(seq))



def insert_del(seq, start, end, zero_based=False):
    '''
    Inserts a deletion, given a sequence and an interval.
    '''
    if end < start:
        print 'Error in deletion! End position must be greater than or equals to Start.'
        return seq
    if not zero_based:
        start -= 1
        end   -= 1
    return seq[:start] + seq[end+1:]



def insert_dup(seq, start, end, zero_based=False):
    '''
    Inserts a duplication, given a sequence and an interval.
    '''
    if end < start:
        print 'Error in duplication! End position must be greater than or equals to Start.'
        return seq
    if not zero_based:
        start -= 1
        end   -= 1
    return seq[:end+1] + seq[start:]



def insert_inv(seq, start, end, zero_based=False):
    if end < start:
        print 'Error in inversion! End position must be greater than or equals to Start.'
        return seq
    if not zero_based:
        start -= 1
        end   -= 1
    return seq[:start] + reverse_complement(seq[start:end+1:]) + seq[end+1:]



def insert_trans(large_seq, small_seq, start, zero_based=False):
    if not zero_based:
        start -= 1
    return large_seq[:start] + small_seq + large_seq[start:]



def insert_indel(seq, start, ref, alt, message):
    '''
    Function to read indel code for pair reference (ref) - alternative (alt)
    and insert the corresponding indel in the correct position ('start') on
    sequence ('seq')
    It returns that resulting sequence and the shift caused by the indel in
    the relative positions of the sequence in the right of the indel.
    '''
    # parsing alternative allele
    i = 0
    while alt[:i+1].isdigit() and i < len(alt):
        i+=1
    subs_seq = alt[i:]
    shift = len(subs_seq) - len(ref) # Value with positions shifted to the right (positive),
                                     # to the left (negative). If shift == 0, nothing changes.
                                     # General formula = positions inserted - positions deleted
    context = 50
    if i > 0:
        number = int(alt[:i])
        if number > 0 and number != len(ref):
            print 'Error! Indel codes of ref[%s] and alt[%s] are not matching!' % (ref, alt)
            # sys.exit(1)
    if i == len(alt): # no letters, just numbers (one deletion, no insertion)
        # print 'deletion! [%s] :: [%s] (pos:%d) <%s>' % (seq[start-1:start-1+len(ref)], ref, start, message)
        seq = insert_del(seq, start, start+number-1)
        shift = -number
    elif i == 0: # no number, just letters (deletion + insertion)
        # if seq[start-1:start-1+len(ref)] != ref:
        # print 'i==0! [%s] != [%s] (pos:%d) <%s>' % (seq[start-1:start-1+len(ref)], ref, start, message)
        # print 'Context (%d:%d) -> [%s]\n' % (start-1-context, start-1+len(ref)+context, seq[start-1-context:start-1+len(ref)+context])
        seq = insert_del(seq, start, start+len(ref)-1)
        seq = insert_trans(seq, subs_seq, start)
    elif number == 0: # no deletion, just insertion
        # print 'number==0! [%s] :: [%s] (pos:%d) <%s>' % (seq[start-1:start-1+len(ref)], ref, start, message)
        seq = insert_trans(seq, subs_seq, start+len(ref)-1)
        shift = len(subs_seq)
    else: # numbers and letters (deletion + insertion)
        # if seq[start-1:start-1+len(ref)] != ref:
        # print 'else [ref=%s|alt=%s]! [%s] :: [%s] (pos:%d) <%s>' % (ref, alt, seq[start-1:start-1+len(ref)], ref, start, message)
        # print 'Context (%d:%d) -> [%s]\n' % (start-1-context, start-1+len(ref)+context, seq[start-1-context:start-1+len(ref)+context])
        seq = insert_del(seq, start, start+len(ref)-1)
        seq = insert_trans(seq, subs_seq, start)
    return seq, shift



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
    # Required
    parser.add_argument('--fasta_input', '-i',  required=True, metavar='input.fasta',  type=file, help='Fasta file to be changed with SVs.')
    parser.add_argument('--fasta_output', '-o', required=True, metavar='output.fasta', type=argparse.FileType('w'), help='Fasta file to be created with SVs.')
    parser.add_argument('--bed', dest='sv_bed', required=True, metavar='SVs.bed', type=file, help='BED file with SVs to be inserted.')
    parser.add_argument('--chrom_lens', required=True, type=file, metavar='chrom_lengths_file', dest='chrom_lens_file', help='Text file with chromosome lengths.')
    parser.add_argument('--chrom', required=True, type=str, metavar='chromosome_name', dest='chromosome_name', help='Chromosome.')
    # Optional
    parser.add_argument('--indels', dest='indels_file', required=False, type=file, help='Text file with indels to be inserted.')
    parser.add_argument('--fasta_label', required=False, type=str, metavar='fasta_label', dest='fasta_label', help='Name to label fasta sequence.')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--overlap', action='store_true')

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

    all_regions = [] # this list will have both large SVs and indels
    for line in lines:
        elems = line.split()
        chrom, start, end, sv_type = elems[:4]
        if chrom == chrom_name:
            if len(elems) > 4:
                translocation_info = elems[4]
                all_regions.append([chrom, int(start), int(end), sv_type.lower(), translocation_info])
            else:
                all_regions.append([chrom, int(start), int(end), sv_type.lower()])

    if args.indels_file:
        lines = args.indels_file.read().split('\n')
        while lines[-1] == '':
            lines.pop()
        for line in lines:
            elems = line.split()
            chrom, start, ref, alt = elems[:4]
            if chrom == chrom_name:
                all_regions.append([chrom, int(start), int(start)+len(ref)-1, 'indel', ref+'__'+alt])

    new_regions = merge_regions(all_regions)
    sv_regions = [] # only large SVs

    for region in new_regions[::-1]:
        if len(region) == 1:
            if region[0][3] in ['dup', 'del', 'inv', 'unbtr', 'baltr']:
                sv_regions.append(region[0])
            elif region[0][3] == 'indel':
                start = region[0][1]
                ref, alt = region[0][4].split('__')
                chrom_seq, shift = insert_indel(chrom_seq, start, ref, alt, 'indel only')
            else:
                print 'Error! Unexpected region type ->', region[0][3]
        elif region[0][3] in ['dup', 'del', 'inv', 'unbtr', 'baltr']:
            large_sv = region[0]
            if large_sv[3] == 'del':
                sv_regions.append(large_sv)
            elif large_sv[3] in ['dup', 'inv']:
                for indel in region[:0:-1]:
                    if indel[3] == 'indel':
                        start = indel[1]
                        ref, alt = indel[4].split('__')
                        chrom_seq, shift = insert_indel(chrom_seq, start, ref, alt, 'indel inside SV')
                        large_sv[2] += shift
                    else:
                        print 'Ignoring overlapping SV -> [%s|%s:%d-%d]. Only larger SV [%s|%s:%d-%d] will be considered.' % \
                               (indel[3], indel[0], indel[1], indel[2], large_sv[3], large_sv[0], large_sv[1], large_sv[2])
                sv_regions.append(large_sv)
        else:
            print 'Error! Unexpected region configuration ->', region


    sv_regions = sorted(sv_regions, key=itemgetter(1), reverse=True)
    # for sv in sv_regions:
    #     print sv

    '''
    new_sv_regions = [sv_regions[0]]
    # Checking intervals overlap
    if not args.overlap:
        if len(sv_regions) > 1:
            for i in range(1, len(sv_regions)):
                interval_overlap = intervals_overlap(sv_regions[i][1:3], sv_regions[i-1][1:3])
                if interval_overlap > 0: # and :
                    # if 
                    print '\n\tError! SV regions are overlapping!'
                    print '\t%s:%d-%d [%s] <--> %s:%d-%d [%s] \n' % (sv_regions[i-1][0], sv_regions[i-1][1], sv_regions[i-1][2], sv_regions[i-1][3],
                                                                     sv_regions[i][0],   sv_regions[i][1],   sv_regions[i][2],   sv_regions[i][3])
                    sys.exit(1)
                else:
                    new_sv_regions.append(sv_regions[i])
    '''


    print 'Inserting SVs in sequence.'
    for sv_region in sv_regions:
        chrom, start, end, sv_type = sv_region[:4]
        if chrom_name == chrom:
            if args.verbose:
                print 'Inserting %s in %s:%s-%s.' % (sv_type, chrom, start, end)
            if len(sv_region) > 4:
                translocation_info = sv_region[4]
            if sv_type.startswith('del'):
                chrom_seq = insert_del(chrom_seq, start, end)
            elif sv_type.startswith('dup'):
                chrom_seq = insert_dup(chrom_seq, start, end)
            elif sv_type.startswith('inv'):
                chrom_seq = insert_inv(chrom_seq, start, end)
            elif sv_type.startswith('unbtr') or sv_type.startswith('baltr'):
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
                trans_seq = get_subseq_from_fasta(trans_chrom_filename, start_line, start_pos, end_line, end_pos).upper()
                if sv_type.startswith('baltr'):
                    chrom_seq = insert_del(chrom_seq, start, end)
                chrom_seq = insert_trans(chrom_seq, trans_seq, start)


    print 'Writing to output file [%s].' % args.fasta_output.name
    write_fasta(chrom_name, args.fasta_output, chrom_seq, fasta_size_per_line, chrom_label)
    print 'Simulations for', chrom_name, 'done.'
    print 'Total size:', LEN_CHROM
    print 


if __name__ == '__main__':
    main()


