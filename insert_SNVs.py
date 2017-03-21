#!/usr/bin/python2
# insert_SNVs.py
# C: Sep 29, 2015
# M: Jan 28, 2016
# A: Leandro Lima <leandrol@usc.edu>


prog_name = 'insert_SNVs.py'


import sys, argparse
from random import randint
from operator import itemgetter

SNVs_inserted = 0


def insert_SNV(REF, ALT, MAF):
    global SNVs_inserted
    rand_number = randint(1, 1000)
    if rand_number <= round(MAF * 1000) and ALT in ['A', 'C', 'G', 'T']:
        SNVs_inserted += 1
        return ALT
    else:
        return REF


def main():

    parser = argparse.ArgumentParser(description='Get arguments to create single-nucleotide variants (SNVs) in a fasta file.', prog=prog_name)

    one_based = 1 # it could be used in case we have zero based positions

    parser.add_argument('--fasta_input',  required=True, metavar='input.fasta',  type=file, help='Input fasta file to be used as reference to receive SNVs.')
    parser.add_argument('--fasta_output', required=True, metavar='output.fasta', dest='output_fasta_file', type=argparse.FileType('w'), help='Output fasta file with random SNVs, based in frequencies.')
    parser.add_argument('--freq_file',    required=True, type=file, help='Text file SNV frequencies.')
    parser.add_argument('--vcf_output',   required=True, type=argparse.FileType('w'), help='VCF file generated with SNVs inserted.')
    parser.add_argument('--chrom', required=True, type=str, metavar='chrom', dest='chromosome_name_vcf', help='Chromosome.')

    args = parser.parse_args()
    
    lines_fasta = args.fasta_input.read().split('\n')
    while lines_fasta[-1] == '':
        lines_fasta.pop()

    if not lines_fasta[0].startswith('>'):
        print 'Invalid fasta file!!'
        sys.exit(1)

    # Frequencies have to be sorted by position
    frequencies = args.freq_file.read().split('\n')
    while frequencies[-1] == '':
        frequencies.pop()

    chrom_name = lines_fasta[0].replace('>', '')
    fasta_size_per_line = len(lines_fasta[1])

    args.output_fasta_file.write('>' + chrom_name + '\n')
    # args.output_fasta_file.write('>' + chrom_name + '_SNVs\n')

    start = 0
    end = start + fasta_size_per_line - 1
    freqs_line = 0 # for frequencies
    CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
    print 'Looking for chromosome %s in frequencies file [%s].' % (args.chromosome_name_vcf, args.freq_file.name)
    while CHROM != args.chromosome_name_vcf:
        freqs_line += 1
        CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
    

    args.vcf_output.write('##fileformat=VCFv4.1\n')
    args.vcf_output.write('##reference=%s\n' % args.fasta_input.name)
    args.vcf_output.write('##contig=<ID=%s>\n' % args.chromosome_name_vcf)
    args.vcf_output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')


    # Reading fasta lines to compare with SNV positions
    for line_fasta in lines_fasta[1:]:
        new_line = line_fasta
        # If there is still frequencies to read
        if freqs_line < len(frequencies):
            # Look for current SNV in fasta line
            while int(POS) <= end + 1:
                if len(REF)*len(ALT) != 1: # inserting only 1-base SNVs
                    freqs_line += 1
                    try:
                        CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                    except:
                        break
                    if CHROM != args.chromosome_name_vcf: # in case finds next chromosome in frequencies file
                        break
                    continue
                elif CHROM == args.chromosome_name_vcf: # confirming chromosome according to input
                    if new_line[int(POS)-start-one_based].upper() == REF:
                        base = insert_SNV(REF, ALT, float(MAF))
                        if new_line[int(POS)-start-one_based].islower():
                            base = base.lower()
                        new_line = new_line[:int(POS)-start-one_based] + base + new_line[int(POS)-start-one_based+1:]
                        if base.upper() != REF:
                            args.vcf_output.write('%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\n' % (CHROM, POS, rsID, REF, ALT))
                    else:
                        pass
                    if int(POS) >= end + 1:
                        freqs_line += 1
                        if freqs_line < len(frequencies):
                            CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                            if CHROM != args.chromosome_name_vcf:
                                break
                        break
                    freqs_line += 1
                    if freqs_line == len(frequencies):
                        break
                    else:
                        CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                        if CHROM != args.chromosome_name_vcf:
                            break
                else:
                    break

        args.output_fasta_file.write(new_line + '\n')
        start = end + 1
        end = start + fasta_size_per_line - 1

    print 'Simulations for', chrom_name, 'done.'
    print SNVs_inserted, 'SNVs were inserted.'
    # print 'Frequency:', 
    print 
    


if __name__ == '__main__':
    main()


