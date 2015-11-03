# simulate_SV_BED.py
# C: Oct  1, 2015
# M: Nov  2, 2015
# A: Leandro Lima <leandrol@usc.edu>


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
    inv_len_filename   = sys.argv[3] # File with duplication lengths (one per line)
    output_bed_file    = sys.argv[4] # Output file name
    chromosome_name    = sys.argv[5] # Chromosome name
    chrom_lens_file    = sys.argv[6] # Chromosome lengths file
    avoid_regions_file = sys.argv[7] # BED file name

    # Usage: python simulate_SV_BED.py [DEL_lengths_file] [DUP_lengths_file] [DUP_lengths_file] [simulated_SVs.bed] [chrom] [chrom_lens_file] [gaps_file]
    # Example: python simulate_SV_BED.py SV_lengths.txt SV_lengths.txt SV_lengths.txt simulated_SVs.bed chrX reference/chrom_lengths_hg19.txt reference/gaps_hg19.txt
    # Example for all chromosomes:
    # mkdir simulated_SVs; for chrom in {1..22} X Y; do python simulate_SV_BED.py SV_lengths.txt SV_lengths.txt simulated_SVs/simulated_SVs_chr$chrom.bed chr$chrom reference/chrom_lengths_hg19.txt reference/gaps_hg19.txt; done

    if chromosome_name.startswith('chr'):
        chromosome_name = chromosome_name[3:]

    distance_between_SVs = 100000
    SV_lens = []

    lines = open(del_len_filename).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    SV_lens = ['del_'+line for line in lines]

    lines = open(dup_len_filename).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    SV_lens += ['dup_'+line for line in lines]

    lines = open(inv_len_filename).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    SV_lens += ['inv_'+line for line in lines]

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
    start = avoid_end + distance_between_SVs
    while SV_lens != []:
        SV = SV_lens[randint(0, len(SV_lens)-1)]
        sv_type, sv_size = SV.split('_')
        end = start + int(sv_size)
        if end > chrom_len:
            print 'Impossible to finish process of creating SVs. Not enough chromosome length (size of %s is %d) for simulated SVs.' % (chromosome_name, chrom_len)
            print 'Exiting.'
            break
        
        if between(start, end, avoid_start-distance_between_SVs) or between(avoid_start-distance_between_SVs, avoid_end+distance_between_SVs, start):
            start = avoid_end + int(1.1*distance_between_SVs)
            avoid_start = regions_to_avoid[0][0]
            avoid_end   = regions_to_avoid[0][1]
            regions_to_avoid = regions_to_avoid[1:]
        else:
            BED_file.write('%s\t%d\t%d\t%s\n' % (chromosome_name, start, end, sv_type))
            if verbose:
                print 'Created ', chromosome_name, str(start), str(end), sv_type
            SV_lens.remove(SV)
            start = end + distance_between_SVs


if __name__ == '__main__':
    main()
