# create_reads.py
# C: Sep 29, 2015
# M: Dec  8, 2015
# A: Leandro Lima <leandrol@usc.edu>


import os, sys, argparse #, linecache
from random import *
# from operator import itemgetter
from insert_SVs import reverse_complement

prog_name = 'create_reads.py'


fastq_qual = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


# min and max alpha for beta distribution
alphaMAX = 15
alphaMIN = 4
qualMAX  = 42
qualMIN  = 30

# Parameters for Illumina quality scores
# alpha_q = 2
beta_q  = 20
# mean_beta_dist = alpha/float(alpha+beta)


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


"""
length = 70
qualities = [min(int(betavariate(alpha_for_short_read_quality(betaMAX, betaMIN, i, n), beta) * avg_base_qual / mean_beta_dist), avg_base_qual) for i in range(length)]
""" 



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



# Not working yet
def write_bam_read(bam_file, name, seq, quality):
    a = pysam.AlignedSegment()
    a.query_name = name # "read_28833_29006_6945"
    a.query_sequence = seq # "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((0,10), (2,1), (0,25))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    a.query_qualities = pysam.fromQualityString(quality)
    a.tags = (("NM", 1),
              ("RG", "L1"))
    bam_file.write(a)



def generate_pair_from_fragment(fragment, read_length):
    # return [fragment[:read_length], ''.join(complement[base] for base in reversed(fragment[-read_length:]))]
    return [fragment[:read_length], reverse_complement(fragment[-read_length:]))]



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
        nonN_regions.add(elem.strip('N'))
        
    return list(nonN_regions)



def main():

    parser = argparse.ArgumentParser(description='This program inserts structural variants from a BED file into a FASTA file.', prog=prog_name)
    parser.add_argument('--fasta_input', '-i',  required=True, metavar='input.fasta',  type=file, help='Fasta file to be changed with SVs.')
    parser.add_argument('--output_file', '-o', required=True, metavar='output.fq|output.bam', type=str, help='Output file for reads. It must finish with .fq/.fastq or .bam (paired-end option will automatically change file names to output1.fq and output2.fq).')
    parser.add_argument('-pe', action='store_true', dest='paired_end', help='Add option to generate paired-end reads.')
    # parser.add_argument('--filetype', required=True, metavar='bam|fastq', type=file, help='Data type of output: unaligned "bam" or "fastq" file.')
    parser.add_argument('--cov', required=True, type=float, metavar='coverage', dest='coverage', help='Average coverage.')
    parser.add_argument('--read_len', required=True, type=int, metavar='avg_read_len', dest='avg_read_len', help='Average read length (length is fix for short reads).')
    parser.add_argument('--ins_rate', required=False, type=float, metavar='insertion_rate', dest='insertion_error_rate', help='Insertion error rate for reads.', default=0.12)
    parser.add_argument('--del_rate', required=False, type=float, metavar='deletion_rate', dest='deletion_error_rate', help='Deletion error rate for reads.', default=0.02)
    parser.add_argument('--snp_rate', required=False, type=float, metavar='snp_rate', dest='snp_error_rate', help='SNP error rate for reads.', default=0.01)
    parser.add_argument('--insert_size', required=False, type=int, metavar='insert_size', help='Insert size for short reads.', default=300)
    parser.add_argument('--insert_sd', required=False, type=int, metavar='insert_sd', help='Insert standard deviation for short reads.', default=50)
    parser.add_argument('--alpha', required=False, type=int, metavar='alpha', help='Alpha for beta distribution of read lengths.', default=2)
    parser.add_argument('--beta', required=False, type=int, metavar='beta', help='Beta for beta distribution of read lengths.', default=10)

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    filetype = args.output_file.split('.')[-1]

    # Checking changing function to write output file, depending on type
    if filetype in ['fastq', 'fq']:
        write_read = write_fastq_read
    elif filetype == 'bam':
        write_read = write_bam_read
        import pysam
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

        # PacBio
        if long_reads:
            read_lens = [min(int(betavariate(alpha, beta) * avg_read_len / mean_beta_dist), len(region)) for i in range(n_reads)]

            # sampling fragments without errors from reference
            reads = generate_reads_from_seq(region.upper().replace('N', ''), read_lens)

            # adding errors to reads and writing it to output file
            for read in reads:
                errors = generate_errors(len(read), args.insertion_error_rate, args.deletion_error_rate, args.snp_error_rate)
                read_with_errors = insert_errors_in_seq(read, errors)
                qualities = create_quality_for_pacbio_read(len(read_with_errors), qual_mean, qual_sd)
                write_read(output1, label + '-' + str(read_id), read_with_errors, qualities)
                read_id += 1

        # short reads
        else:
            # individual reads
            if paired_end:

                # insert + pair
                read_lens = [min(2*avg_read_len + int(gauss(args.insert_size, args.insert_sd)), len(region)) for i in range(n_reads)]

                # sampling fragments without errors from reference
                fragments = generate_reads_from_seq(region.upper().replace('N', ''), read_lens)

                # adding errors to reads and writing it to output file
                for fragment in fragments:
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
                reads = generate_reads_from_seq(region.upper().replace('N', ''), read_lens)

                # adding errors to reads and writing it to output file
                for read in reads:
                    errors = generate_errors(len(read), args.insertion_error_rate, args.deletion_error_rate, args.snp_error_rate)
                    read_with_errors = insert_errors_in_seq(read, errors)
                    qualities = create_quality_for_illumina_read(len(read_with_errors))
                    write_read(output1, label + '-' + str(read_id), read_with_errors, qualities)
                    read_id += 1



    output1.close()
    if paired_end:
        output2.close()



if __name__ == '__main__':
    main()


"""

header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                   {'LN': 1584, 'SN': 'chr2'}] }

tmpfilename = 'test.bam'
outfile = pysam.AlignmentFile(tmpfilename, "wb", header=header)
a = pysam.AlignedSegment()
a.query_name = "read_28833_29006_6945"
a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
a.flag = 99
a.reference_id = 0
a.reference_start = 32
a.mapping_quality = 20
a.cigar = ((0,10), (2,1), (0,25))
a.next_reference_id = 0
a.next_reference_start=199
a.template_length=167
a.query_qualities = pysam.fromQualityString("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
a.tags = (("NM", 1),
          ("RG", "L1"))
outfile.write(a)
outfile.write(a)
outfile.close()
"""

