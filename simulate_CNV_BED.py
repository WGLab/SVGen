# simulate_CNV_BED.py
# C: Oct  1, 2015
# M: Oct 19, 2015
# A: Leandro Lima <leandrol@usc.edu> / <llima@ime.usp.br>


import sys
from random import randint
from operator import itemgetter



def between(start, end, value):
    if (start <= value and value <= end) or (end <= value and value <= start):
        return True
    else:
        return False


def main():

    verbose = True # make it an option

    del_len_filename   = sys.argv[1] # File with deletion lengths (one per line)
    dup_len_filename   = sys.argv[2] # File with duplication lengths (one per line)
    output_bed_file    = sys.argv[3] # Output file name
    chromosome_name    = sys.argv[4] # Chromosome name
    chrom_lens_file    = sys.argv[5] # Chromosome lengths file
    avoid_regions_file = sys.argv[6] # BED file name
    # genome_version     = sys.argv[7] # genome version (hg18, hg19, hg38)

    # Usage: python simulate_CNV_BED.py [DEL_lengths_file] [DUP_lengths_file] [simulated_CNVs.bed] [chrom] [chrom_lens_file] [gaps_file]
    # Example: python simulate_CNV_BED.py CNV_lengths.txt CNV_lengths.txt simulated_CNVs.bed chrX reference/chrom_lengths_hg19.txt reference/gaps_hg19.txt

    if chromosome_name.startswith('chr'):
        chromosome_name = chromosome_name[3:]

    distance_between_CNVs = 1000000
    CNV_lens = []

    lines = open(del_len_filename).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    CNV_lens = ['del_'+line for line in lines]

    lines = open(dup_len_filename).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    CNV_lens += ['dup_'+line for line in lines]

    # Avoiding gap regions (centromeres, telomeres, etc.): https://www.biostars.org/p/2349/#2351

    lines = open(avoid_regions_file).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    while lines[0].startswith('#'):
        lines = lines[1:]

    regions_to_avoid = []
    for line in lines:
        chrom, start, end = line.split()[1:4]
        if chrom.replace('chr', '') == chromosome_name:
            regions_to_avoid.append([int(start), int(end)])

    regions_to_avoid = sorted(regions_to_avoid, key=itemgetter(0), reverse=False)
    avoid_start = regions_to_avoid[0][0]
    avoid_end   = regions_to_avoid[0][1]
    regions_to_avoid = regions_to_avoid[1:]

    chrom_lens = {}
    with open(chrom_lens_file) as chrom_lens_file_lines:
        for line in chrom_lens_file_lines:
            try:
                chrom, chrom_len = line.split()
                chrom_lens[chrom] = int(chrom_len)
            except:
                pass


    chrom_len = chrom_lens[chromosome_name]

    BED_file = open(output_bed_file, 'w')
    start = avoid_end + distance_between_CNVs
    while CNV_lens != []:
        CNV = CNV_lens[randint(0, len(CNV_lens)-1)]
        cnv_type, cnv_size = CNV.split('_')
        end = start + int(cnv_size)
        if end > chrom_len:
            print 'Impossible to finish process of creating CNVs. Not enough chromosome length for simulated CNVs.'
        
        if between(start, end, avoid_start-distance_between_CNVs) or between(avoid_start-distance_between_CNVs, avoid_end+distance_between_CNVs, start):
            start = avoid_end + int(1.1*distance_between_CNVs)
            avoid_start = regions_to_avoid[0][0]
            avoid_end   = regions_to_avoid[0][1]
            regions_to_avoid = regions_to_avoid[1:]
        else:
            BED_file.write('%s\t%d\t%d\t%s\n' % (chromosome_name, start, end, cnv_type))
            if verbose:
                print 'created ', chromosome_name, str(start), str(end), cnv_type
            CNV_lens.remove(CNV)
            start = end + distance_between_CNVs


if __name__ == '__main__':
    main()
