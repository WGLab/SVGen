#!/usr/bin/python2
# create_reads.py
# C: Sep 29, 2015
# M: May  3, 2016
# A: Leandro Lima <leandrol@usc.edu>


import os, sys, argparse, cPickle #, linecache
from random import *
# from operator import itemgetter
from insert_SVs import reverse_complement
import numpy as np

prog_name = 'create_reads.py'


fastq_qual = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# Creating GC->probability dictionary
gc = [i/100. for i in range(101)]
cov  = [(((62.9390 + 390.7907 * i + (-439.2405) * i**2) - 144.80) / 38.86921) for i in gc]
cov = [i-min(cov) for i in cov]
cov = [i/max(cov) for i in cov]
gc_cov = {}
gc_cov.update(zip(gc, cov))

# min and max alpha for beta distribution
alphaMAX = 15
alphaMIN = 4

# max and min random data of quality score
qualMAX  = 45
qualMIN  = 25

# Parameters for Illumina quality scores
# alpha_q = 2
beta_q  = 20

# Bin szie for GC content
gc_region_size = 1000



def gc_content(sequence):
    sequence = sequence.upper()
    A_count = sequence.count('A')
    C_count = sequence.count('C')
    G_count = sequence.count('G')
    T_count = sequence.count('T')
    gc = float(G_count + C_count)/(A_count+C_count+G_count+T_count)
    return gc



def gc_intervals(sequence, size, min_bases):
    intervals = [[start, start+size-1, gc_content(sequence[start-1:start+size-1])] for start in range(1, len(sequence)+1, size)
                 if len(sequence[start-1:start+size-1])-sequence[start-1:start+size-1].count('N') >= min_bases]
    return intervals



def prob_intervals(sequence, size, min_bases):
    """
    The same as the function 'gc_intervals', but instead of returning
    GC-content, it directly returns the probabilities
    """
    intervals = [[start, start+size-1, gc_cov[round(gc_content(sequence[start-1:start+size-1]),2)]] for start in range(1, len(sequence)+1, size)
                 if len(sequence[start-1:start+size-1])-sequence[start-1:start+size-1].count('N') >= min_bases]
    return intervals



def gc_to_prob(intervals):
    SUM = sum([interval[2] for interval in intervals])
    interval_probs = [[interval[0], interval[1], interval[2]/SUM] for interval in intervals]
    return interval_probs



def alpha_for_short_read_quality(alphaMAX, alphaMIN, i, n):
    return (alphaMIN-alphaMAX)*(i-1)/(n-1) + alphaMAX



def create_quality_for_pacbio_read(length, mean, sd):
    qualities = [fastq_qual[min(max(0, int(gauss(mean, sd))), 93)] for i in range(length)]
    return ''.join(qualities)



def quality_by_beta_dist(beta_output, qualMIN, qualMAX):
    return qualMAX - (1-beta_output)*(qualMAX-qualMIN)



def create_quality_for_illumina_read(length):
    beta_outputs = [betavariate(alpha_for_short_read_quality(alphaMAX, alphaMIN, i, length), beta_q) for i in range(length)]
    qualities    = [fastq_qual[int(quality_by_beta_dist(beta_output, qualMIN, qualMAX))] for beta_output in beta_outputs]
    return ''.join(qualities)



def generate_snp(base):
    if base == 'A':
        return choice(['C', 'G', 'T'])
    if base == 'C':
        return choice(['A', 'G', 'T'])
    if base == 'G':
        return choice(['A', 'C', 'T'])
    if base == 'T':
        return choice(['A', 'C', 'G'])



def return_bases_by_error(base, error):
    if error == 'i': # insertion
        return base + choice(['A', 'C', 'G', 'T'])
    if error == 'd': # deletion
        return ''
    if error == 's': # SNP
        return generate_snp(base)
    


def generate_errors(length, insertion_rate, deletion_rate, snp_rate):
    '''
    Given a sequence length and common errors rates, this functions returns
    a list of errors with the following code:
    '' : no error
    'i': insertion
    'd': deletion
    's': SNP
    '''
    number_of_insertions = int(length * insertion_rate)
    number_of_deletions  = int(length * deletion_rate)
    number_of_snps       = int(length * snp_rate)

    # adding more randomness in case rates are too low
    if number_of_insertions == 0 and not insertion_rate == 0:
        if randint(1, 10000) <= insertion_rate * 10000:
            number_of_insertions = 1

    if number_of_deletions == 0 and not deletion_rate == 0:
        if randint(1, 10000) <= deletion_rate * 10000:
            number_of_deletions = 1

    if number_of_snps == 0 and not snp_rate == 0:
        if randint(1, 10000) <= snp_rate * 10000:
            number_of_snps = 1

    number_of_errors = number_of_insertions + number_of_deletions + number_of_snps

    error_positions = sample(range(length), number_of_errors)
    insertion_positions = error_positions[:number_of_insertions]
    deletion_positions  = error_positions[number_of_insertions:number_of_insertions+number_of_deletions]
    snp_positions       = error_positions[number_of_insertions+number_of_deletions:]

    errors = length * ['']

    position = 0
    i = 0
    while i < number_of_insertions:
        errors[error_positions[position]] = 'i'
        i+=1
        position+=1

    i=0
    while i < number_of_deletions:
        errors[error_positions[position]] = 'd'
        i+=1
        position+=1

    i=0
    while i < number_of_snps:
        errors[error_positions[position]] = 's'
        i+=1
        position+=1

    return errors



def insert_errors_in_seq(seq, errors):
    '''
    This function receives a sequence (string) and a list of errors.
    The output is another sequence, with the corresponding errors in the list.
    '''
    return ''.join([seq[i] if errors[i] == '' else return_bases_by_error(seq[i], errors[i]) for i in range(len(seq))])



def write_fastq_read(fastq_file, name, seq, quality):
    '''
    This function receives a fastq file and a read, with name, sequence and quality scores.
    It writes the read (4 lines) to the corresponding file.
    '''
    fastq_file.write('@%s\n' % name)
    fastq_file.write('%s\n'  % seq)
    fastq_file.write('+\n')
    fastq_file.write('%s\n' % quality)



def generate_pair_from_fragment(fragment, read_length):
    # return [fragment[:read_length], ''.join(complement[base] for base in reversed(fragment[-read_length:]))]
    return [fragment[:read_length], reverse_complement(fragment[-read_length:])]

"""
def calculate_coverage(regions, gc):
    for i in range(len(regions)):

    ((((62.9390 + 390.7907 * gc + (-439.2405) * gc ^ 2) - 144.80) / 38.86921) + 10)
"""


def generate_reads_from_seq_by_gc(seq, read_lens):
    '''
    This function gets a sequence (string) and a list of read lengths
    and returns a list of reads distributed according to GC content
    throughout the sequence with the corresponding read lengths. Some 
    reads may be clipped if located in the end of the sequence.
    '''
    n_reads = len(read_lens)
    exp_start_distance = len(seq) / len(read_lens)
    min_nonN_bases = 50
    intervals = prob_intervals(seq, gc_region_size, min_nonN_bases)
    intervals = gc_to_prob(intervals)
    # print intervals[:100]
    region_probs = [interval[2] for interval in intervals]
    # regions = [interval[:2] for interval in intervals]
    index_intervals = np.random.choice(range(len(intervals)), n_reads, p=region_probs) # np.random.choice(4, 1000, p=[0.1, 0.2, 0.3, 0.4])

    starts = [randint(intervals[i][0], intervals[i][1]) for i in index_intervals]
    # starts = range(0, len(seq), exp_start_distance)
    ends = [starts[i]+read_lens[i] for i in range(n_reads)]
    reads = []
    for i in range(n_reads):
        if ends[i] > len(seq):
            reads.append(seq[starts[i]:])
        else:
            reads.append(seq[starts[i]:ends[i]])
    return reads


def generate_reads_from_seq(seq, read_lens):
    '''
    This function gets a sequence (string) and a list of read lengths
    and returns a list of reads distributed equally throughout the 
    sequence with the corresponding read lengths. Some reads may be
    clipped if located in the end of the sequence.
    '''
    n_reads = len(read_lens)
    exp_start_distance = len(seq) / len(read_lens)
    starts = range(0, len(seq), exp_start_distance)
    ends = [starts[i]+read_lens[i] for i in range(n_reads)]
    reads = []
    for i in range(n_reads):
        if ends[i] > len(seq):
            reads.append(seq[starts[i]:])
        else:
            reads.append(seq[starts[i]:ends[i]])
    return reads



def get_nonN_regions(seq):
    '''
    Given a sequence with multiple 'N' regions, this function
    returns a list of regions without multiple N's (short N 
    regions could be present in output, for not being large gaps)
    '''
    # lines = open('mother_hg19_chr2.fa').read().split('\n')
    # seq = ''.join(lines[1:])
    nonN_regions = set(seq.split('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'))
    nonN_regions.remove('')
    
    for elem in nonN_regions:
        nonN_regions.remove(elem)
        # print elem.strip('N')[:10] + ' ... ' + elem.strip('N')[-10:] + ' - total: ' + str(len(elem.strip('N')))
        # if len(elem.strip('N')) < 10:
        #     print elem
        nonN_regions.add(elem.strip('N'))
        
    return list(nonN_regions)



def main():

    parser = argparse.ArgumentParser(description='This program inserts structural variants from a BED file into a FASTA file.', prog=prog_name)
    # Required
    parser.add_argument('--fasta_input', '-i',  required=True, metavar='input.fasta',  type=file, help='Fasta file to be changed with SVs.')
    parser.add_argument('--output_file', '-o', required=True, metavar='output.fq|output.bam', type=str, help='Output file for reads. It must finish with .fq/.fastq or .bam (paired-end option will automatically change file names to output1.fq and output2.fq).')
    parser.add_argument('--cov', required=True, type=float, metavar='coverage', dest='coverage', help='Average coverage.')
    parser.add_argument('--read_len', required=True, type=int, metavar='avg_read_len', dest='avg_read_len', help='Average read length (length is fix for short reads).')
    # Optional
    parser.add_argument('-pe', action='store_true', dest='paired_end', help='Add option to generate paired-end reads.')
    parser.add_argument('--ins_rate', required=False, type=float, metavar='insertion_rate', dest='insertion_error_rate', help='Insertion error rate for reads.', default=0.12)
    parser.add_argument('--del_rate', required=False, type=float, metavar='deletion_rate', dest='deletion_error_rate', help='Deletion error rate for reads.', default=0.02)
    parser.add_argument('--snp_rate', required=False, type=float, metavar='snp_rate', dest='snp_error_rate', help='SNP error rate for reads.', default=0.01)
    parser.add_argument('--insert_size', required=False, type=int, metavar='insert_size', help='Insert size for short reads.', default=300)
    parser.add_argument('--insert_sd', required=False, type=int, metavar='insert_sd', help='Insert standard deviation for short reads.', default=50)
    parser.add_argument('--alpha', required=False, type=int, metavar='alpha', help='Alpha for beta distribution of read lengths.', default=2)
    parser.add_argument('--beta', required=False, type=int, metavar='beta', help='Beta for beta distribution of read lengths.', default=10)
    parser.add_argument('--read_label', required=False, type=str, help='Label to add in each read.')
    parser.add_argument('--fast_sample', required=False, type=int, metavar='fast_sample', default=1000,
                        help='Number of quality score strings, or sets of errors, for reads. Setting a number for this variable will make the process of creating quality scores faster.')

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    filetype = args.output_file.split('.')[-1]

    # Checking changing function to write output file, depending on type
    if filetype in ['fastq', 'fq']:
        write_read = write_fastq_read
    # elif filetype == 'bam':
        # write_read = write_bam_read
        # import pysam
    else:
        print '\nError: output file must be a .fq/.fastq or .bam file.'
        sys.exit(1)

    paired_end = args.paired_end

    if paired_end:
        prefix = '.'.join(args.output_file.split('.')[:-1])
        output1 = open(prefix + '1.' + filetype, 'w')
        output2 = open(prefix + '2.' + filetype, 'w')
    else:
        output1 = open(args.output_file, 'w')

    lines = args.fasta_input.read().split('\n')
    label = lines[0]
    if not label.startswith('>'):
        print 'Invalid fasta file!!'
        sys.exit(1)
    else:
        if not args.read_label is None:
            label = args.read_label + '_' + label[1:].split()[0]
        else:
            label = label[1:].split()[0]

    if args.verbose:
        print 'Finding non-N regions.'

    seq = ''.join(lines[1:])
    nonN_regions = get_nonN_regions(seq)

    if args.verbose:
        print len(nonN_regions), 'non-N regions were found.'

    alpha = args.alpha
    beta  = args.beta
    avg_read_len = args.avg_read_len
    mean_beta_dist = alpha/float(alpha+beta)


    if avg_read_len > 500:
        long_reads = True
        qual_mean = 30
        qual_sd = 10
    else:
        long_reads = False

    if not args.fast_sample is None:
        if long_reads:
            read_qualities = [create_quality_for_pacbio_read(avg_read_len, qual_mean, qual_sd) for i in range(int(args.fast_sample))]
        else:
            read_qualities = [create_quality_for_illumina_read(avg_read_len) for i in range(int(args.fast_sample))]
        reads_errors = [generate_errors(avg_read_len, args.insertion_error_rate, args.deletion_error_rate, args.snp_error_rate) for i in range(int(args.fast_sample))]
        # print read_qualities

    if args.verbose:
        print 'Creating reads.'

    region_count = 0
    read_id = 1

    for region in nonN_regions:
        if args.verbose:
            print 'Region %3d of %3d (length = %8s).' % (region_count+1, len(nonN_regions), len(region))
            region_count += 1

        # region = ''.join([choice(['A', 'C', 'G', 'T']) for i in range(1000)])
        n_reads = int(args.coverage * len(region) / avg_read_len)
        if paired_end:
            n_reads /= 2


        # PacBio
        if long_reads:
            if n_reads > 0:
                read_lens = [min(int(betavariate(alpha, beta) * avg_read_len / mean_beta_dist), len(region)) for i in range(n_reads)]

                # sampling fragments without errors from reference
                reads = generate_reads_from_seq_by_gc(region.upper().replace('N', ''), read_lens)

                # adding errors to reads and writing it to output file
                for read in reads:
                    if args.fast_sample is None:
                        qualities = create_quality_for_pacbio_read(len(read_with_errors), qual_mean, qual_sd)
                        errors = generate_errors(len(read), args.insertion_error_rate, args.deletion_error_rate, args.snp_error_rate)
                    else:
                        qualities = choice(read_qualities)
                        errors = choice(reads_errors)
                    read_with_errors = insert_errors_in_seq(read, errors)
                    write_read(output1, label + '-' + str(read_id), read_with_errors, qualities)
                    read_id += 1

        # short reads
        else:
            # individual reads
            if paired_end:

                # insert + pair
                read_lens = [min(2*avg_read_len + int(gauss(args.insert_size, args.insert_sd)), len(region)) for i in range(n_reads)]

                # sampling fragments without errors from reference
                fragments = generate_reads_from_seq_by_gc(region.upper().replace('N', ''), read_lens)

                # adding errors to reads and writing it to output file
                for fragment in fragments:
                    if len(fragment) > 2.5*avg_read_len:
                        read1, read2 = generate_pair_from_fragment(fragment, avg_read_len)
                        errors1 = generate_errors(len(read1), args.insertion_error_rate, args.deletion_error_rate, args.snp_error_rate)
                        errors2 = generate_errors(len(read2), args.insertion_error_rate, args.deletion_error_rate, args.snp_error_rate)
                        read_with_errors1 = insert_errors_in_seq(read1, errors1)
                        read_with_errors2 = insert_errors_in_seq(read2, errors2)
                        qualities1 = create_quality_for_illumina_read(len(read_with_errors1))
                        qualities2 = create_quality_for_illumina_read(len(read_with_errors2))
                        write_read(output1, label + '-' + str(read_id) + '/1', read_with_errors1, qualities1)
                        write_read(output2, label + '-' + str(read_id) + '/2', read_with_errors2, qualities2)
                        read_id += 1


            else:
                read_lens = [avg_read_len] * n_reads

                # sampling fragments without errors from reference
                reads = generate_reads_from_seq_by_gc(region.upper().replace('N', ''), read_lens)

                # adding errors to reads and writing it to output file
                for read in reads:
                    errors = generate_errors(len(read), args.insertion_error_rate, args.deletion_error_rate, args.snp_error_rate)
                    read_with_errors = insert_errors_in_seq(read, errors)
                    qualities = create_quality_for_illumina_read(len(read_with_errors))
                    write_read(output1, label + '-' + str(read_id), read_with_errors, qualities)
                    read_id += 1

    # closing files
    output1.close()
    if paired_end:
        output2.close()


if __name__ == '__main__':
    main()

