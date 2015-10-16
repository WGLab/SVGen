# simulate_CNV_BED.py
# C: Oct  1, 2015
# M: Oct 16, 2015
# A: Leandro Lima <llima@ime.usp.br>


import sys
from random import randint
from operator import itemgetter


chrom_lens = {

    "hg19": {
        "chr1"  : 249250621,
        "chr2"  : 243199373,
        "chr3"  : 198022430,
        "chr4"  : 191154276,
        "chr5"  : 180915260,
        "chr6"  : 171115067,
        "chr7"  : 159138663,
        "chr8"  : 146364022,
        "chr9"  : 141213431,
        "chr10" : 135534747,
        "chr11" : 135006516,
        "chr12" : 133851895,
        "chr13" : 115169878,
        "chr14" : 107349540,
        "chr15" : 102531392,
        "chr16" : 90354753,
        "chr17" : 81195210,
        "chr18" : 78077248,
        "chr19" : 59128983,
        "chr20" : 63025520,
        "chr21" : 48129895,
        "chr22" : 51304566,
    }
}

def between(start, end, value):
    if (start <= value and value <= end) or (end <= value and value <= start):
        return True
    else:
        return False


def main():

    del_len_filename   = sys.argv[1] # File with deletion lengths (one per line)
    dup_len_filename   = sys.argv[2] # File with duplication lengths (one per line)
    output_bed_file    = sys.argv[3] # Output file name
    chromosome_name    = sys.argv[4] # Chromosome name
    avoid_regions_file = sys.argv[5] # BED file name

    if not chromosome_name.startswith('chr'):
        chromosome_name = 'chr' + chromosome_name
    
    # Change it later to make it flexible
    genome_version = 'hg19'

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
        if chrom == chromosome_name:
            regions_to_avoid.append([int(start), int(end)])

    regions_to_avoid = sorted(regions_to_avoid, key=itemgetter(0), reverse=False)
    avoid_start = regions_to_avoid[0][0]
    avoid_end   = regions_to_avoid[0][1]
    regions_to_avoid = regions_to_avoid[1:]

    chrom_len = chrom_lens[genome_version][chromosome_name]

    BED_file = open(output_bed_file, 'w')
    start = avoid_end + distance_between_CNVs
    while CNV_lens != []:
        CNV = CNV_lens[randint(0, len(CNV_lens)-1)]
        cnv_type, cnv_size = CNV.split('_')
        end = start + int(cnv_size)
        if end > chrom_lens[genome_version][chromosome_name]:
            print 'Impossible to finish process of creating CNVs. Not enough chromosome length for simulated CNVs.'
        
        if between(start, end, avoid_start-distance_between_CNVs) or between(avoid_start-distance_between_CNVs, avoid_end+distance_between_CNVs, start):
            start = avoid_end + int(1.1*distance_between_CNVs)
            avoid_start = regions_to_avoid[0][0]
            avoid_end   = regions_to_avoid[0][1]
            regions_to_avoid = regions_to_avoid[1:]
        else:
            BED_file.write('%s\t%d\t%d\t%s\n' % (chromosome_name, start, end, cnv_type))
            #print 'created ', chromosome_name, str(start), str(end), cnv_type
            CNV_lens.remove(CNV)
            start = end + distance_between_CNVs


if __name__ == '__main__':
    main()
