# simulate_SNVs.py
# C: Sep 29, 2015
# M: Oct 19, 2015
# A: Leandro Lima <leandrol@usc.edu>

############################################################
###
### .. notes ..
### make it parallel (at least the SNV portion)
###
############################################################

import sys
from random import randint
from operator import itemgetter

SNVs_inserted = 0


def insert_SNV(REF, ALT, MAF):
    global SNVs_inserted
    rand_number = randint(1, 1000)
    if rand_number <= round(MAF * 1000):
        SNVs_inserted += 1
        # print 'SNVs_inserted:', SNVs_inserted
        return ALT
    else:
        return REF



def main():

    one_based = 1
    
    fasta_filename      = sys.argv[1]
    freq_filename       = sys.argv[2]
    output_vcf_filename = sys.argv[3]

    output_vcf = open(output_vcf_filename, 'w')

    lines_fasta = open(fasta_filename).read().split('\n')
    while lines_fasta[-1] == '':
        lines_fasta.pop()

    if not lines_fasta[0].startswith('>'):
        print 'Invalid fasta file!!'
        sys.exit(1)

    # /home/llima/references/hg19_1000g2015aug/hg19_ALL.sites.2015_08.chr*.txt
    # Downloaded from http://www.openbioinformatics.org/annovar/download/hg19_1000g2015aug.zip
    # Make sure it's sorted by position
    frequencies = open(freq_filename).read().split('\n')
    while frequencies[-1] == '':
        frequencies.pop()


    chrom_name = lines_fasta[0].replace('>', '')
    fasta_size_per_line = len(lines_fasta[1])

    output_fasta_file = open(chrom_name + '_SNVs.fa', 'w')
    output_fasta_file.write('>' + chrom_name + '_SNVs\n')

    start = 0
    end = start + fasta_size_per_line - 1
    freqs_line = 0 # for frequencies
    CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split() # The program is assuming that the chromosome is the same. Fix this later!
    print CHROM, POS, REF, ALT, MAF, rsID

    for line_fasta in lines_fasta[1:]:
        # print 'start', start, 'end', end
        new_line = line_fasta
        if freqs_line < len(frequencies):
            while int(POS) <= end + 1 and len(REF)*len(ALT) == 1:
                if new_line[int(POS)-start-one_based].upper() == REF:
                    base = insert_SNV(REF, ALT, float(MAF))
                    if new_line[int(POS)-start-one_based].islower():
                        base = base.lower()
                    new_line = new_line[:int(POS)-start-one_based] + base + new_line[int(POS)-start-one_based+1:]
                    if base.upper() != REF:
                        output_vcf.write('%s\t%s\t%s\t%s\t.\t%s\n' % (CHROM, POS, REF, ALT, rsID))
                    # print 'new_line', new_line, ' = ', new_line[:int(POS)-start-one_based], base, new_line[int(POS)-start-one_based+1:], '(', line_fasta, ')'
                # if int(POS) == end + 1:
                if int(POS) >= end + 1:
                    freqs_line += 1
                    if freqs_line < len(frequencies):
                        CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                        print CHROM, POS, REF, ALT, MAF, rsID
                    print 'breaking... (pos >= end)'
                    break
                freqs_line += 1
                if freqs_line == len(frequencies):
                    break
                else:
                    CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                    # print CHROM, POS, REF, ALT, MAF, rsID

        output_fasta_file.write(new_line + '\n')
        start = end + 1
        end = start + fasta_size_per_line - 1

    print 'Simulations for', chrom_name, 'done.'
    print SNVs_inserted, 'SNVs were inserted.'
    print 
    


if __name__ == '__main__':
    main()


