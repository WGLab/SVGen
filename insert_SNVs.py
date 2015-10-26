# insert_SNVs.py
# C: Sep 29, 2015
# M: Oct 26, 2015
# A: Leandro Lima <leandrol@usc.edu>


import sys
from random import randint
from operator import itemgetter

SNVs_inserted = 0


def insert_SNV(REF, ALT, MAF):
    global SNVs_inserted
    rand_number = randint(1, 1000)
    if rand_number <= round(MAF * 1000):
        SNVs_inserted += 1
        return ALT
    else:
        return REF


def main():

    one_based = 1 # it could be used in case we have zero based positions
    
    fasta_filename      = sys.argv[1]
    freq_filename       = sys.argv[2]
    chromosome_name_vcf = sys.argv[3]
    output_vcf_filename = sys.argv[4]

    output_vcf = open(output_vcf_filename, 'w')

    lines_fasta = open(fasta_filename).read().split('\n')
    while lines_fasta[-1] == '':
        lines_fasta.pop()

    if not lines_fasta[0].startswith('>'):
        print 'Invalid fasta file!!'
        sys.exit(1)

    # Frequencies have to be sorted by position
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
    CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
    print 'Looking for chromosome in frequencies file.'
    while CHROM != chromosome_name_vcf:
        freqs_line += 1
        CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
        
    # print CHROM, POS, REF, ALT, MAF, rsID


    # Reading fasta lines to compare with SNV positions
    for line_fasta in lines_fasta[1:]:
        new_line = line_fasta
        # If there is still frequencies to read
        if freqs_line < len(frequencies):
            # Look for current SNV in fasta line
            while int(POS) <= end + 1:
                if len(REF)*len(ALT) != 1: # inserting only 1-base SNVs
                    freqs_line += 1
                    CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                    if CHROM != chromosome_name_vcf: # in case finds next chromosome in frequencies file
                        break
                    continue
                elif CHROM == chromosome_name_vcf: # confirming chromosome according to input
                    if new_line[int(POS)-start-one_based].upper() == REF:
                        base = insert_SNV(REF, ALT, float(MAF))
                        if new_line[int(POS)-start-one_based].islower():
                            base = base.lower()
                        new_line = new_line[:int(POS)-start-one_based] + base + new_line[int(POS)-start-one_based+1:]
                        if base.upper() != REF:
                            output_vcf.write('%s\t%s\t%s\t%s\t.\t%s\n' % (CHROM, POS, REF, ALT, rsID))
                    else:
                        # Print in case fasta position is different of freq. position
                        print 'There is something wrong with position', POS, ' - fasta =', new_line[int(POS)-start-one_based], ' - vcf =', REF
                    if int(POS) >= end + 1:
                        freqs_line += 1
                        if freqs_line < len(frequencies):
                            CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                            if CHROM != chromosome_name_vcf:
                                break
                        break
                    freqs_line += 1
                    if freqs_line == len(frequencies):
                        break
                    else:
                        CHROM, POS, REF, ALT, MAF, rsID = frequencies[freqs_line].split()
                        # print CHROM, POS, REF, ALT, MAF, rsID
                        if CHROM != chromosome_name_vcf:
                            break
                else:
                    break

        output_fasta_file.write(new_line + '\n')
        start = end + 1
        end = start + fasta_size_per_line - 1

    print 'Simulations for', chrom_name, 'done.'
    print SNVs_inserted, 'SNVs were inserted.'
    print 
    


if __name__ == '__main__':
    main()


