# simulate_SV_BED.py
# C: Oct  1, 2015
# M: Nov  6, 2015
# A: Leandro Lima <leandrol@usc.edu>


import sys, argparse
from random import randint
from operator import itemgetter

prog_name = 'simulate_SV_BED.py'


def between(start, end, value):
    if (start <= value and value <= end) or (end <= value and value <= start):
        return True
    else:
        return False


def main():

    parser = argparse.ArgumentParser(description='Get arguments to create random structural variant regions in a BED file.', prog=prog_name)
    parser.add_argument('--del_lens',   nargs='?', type=file, metavar='del_lengths_file',   dest='del_len_filename',   help='Text file with deletion lengths.')
    parser.add_argument('--dup_lens',   nargs='?', type=file, metavar='dup_lengths_file',   dest='dup_len_filename',   help='Text file with duplication lengths.')
    parser.add_argument('--inv_lens',   nargs='?', type=file, metavar='inv_lengths_file',   dest='inv_len_filename',   help='Text file with inversion lengths.')
    parser.add_argument('--trans_lens', nargs='?', type=file, metavar='trans_lengths_file', dest='trans_len_filename', help='Text file with translocation lengths.')
    
    parser.add_argument('--chrom_lens', required=True, type=file, metavar='chrom_lengths_file', dest='chrom_lens_file', help='Text file with chromosome lengths.')
    parser.add_argument('--output', '-o', required=True, type=argparse.FileType('w'), metavar='output_bed_file', dest='output_bed_file', help='BED output file.')
    parser.add_argument('--gaps', required=True, type=file, metavar='gaps_file', dest='avoid_regions_file', help='BED file with regions to avoid (centromeres and telomeres).')
    parser.add_argument('--chrom', required=True, type=str, metavar='chromosome_name', dest='chromosome_name', help='Chromosome.')
    parser.add_argument('--distance', '-d', type=int, metavar='distance_between_SVs', dest='distance_between_SVs', help='Distance between SVs.', default=100000)

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    # Checking if at least one type of variant was provided
    if args.del_len_filename is None and args.dup_len_filename is None and args.inv_len_filename is None and args.trans_len_filename is None:
        print '\n\n\tError! At least one file with variant (deletion, duplication, inversion or translocation) lengths has to be provided.\n\n'
        sys.exit(1)


    # Example:
    #       python simulate_SV_BED.py --del_lens SV_lengths.txt --dup_lens SV_lengths.txt -o simulated_SVs_chr7.bed --chrom chr7 --chrom_lens reference/chrom_lengths_hg19.txt --gaps reference/gaps_hg19.txt
    #
    # Example for all chromosomes:
    #       mkdir simulated_SVs;
    #       for chrom in {1..22} X Y; do python simulate_SV_BED.py --del_lens SV_lengths.txt --dup_lens SV_lengths.txt -o simulated_SVs/simulated_SVs_chr$chrom.bed --chrom chr$chrom --chrom_lens reference/chrom_lengths_hg19.txt --gaps reference/gaps_hg19.txt; done

    if args.chromosome_name.startswith('chr'):
        args.chromosome_name = args.chromosome_name[3:]

    SV_lens = []

    if not args.del_len_filename is None:
        lines = args.del_len_filename.read().split('\n')
        while lines[-1] == '':
            lines.pop()
        SV_lens = ['del_'+line for line in lines]

    if not args.dup_len_filename is None:
        lines = args.dup_len_filename.read().split('\n')
        while lines[-1] == '':
            lines.pop()
        SV_lens += ['dup_'+line for line in lines]

    if not args.inv_len_filename is None:
        lines = args.inv_len_filename.read().split('\n')
        while lines[-1] == '':
            lines.pop()
        SV_lens += ['inv_'+line for line in lines]

    if not args.trans_len_filename is None:
        lines = args.trans_len_filename.read().split('\n')
        while lines[-1] == '':
            lines.pop()
        SV_lens += ['trans_'+line for line in lines]

    # Avoiding gap regions (centromeres, telomeres, etc.): https://www.biostars.org/p/2349/#2351

    lines = args.avoid_regions_file.read().split('\n')
    while lines[-1] == '':
        lines.pop()

    while lines[0].startswith('#'):
        lines = lines[1:]

    regions_to_avoid = []
    for line in lines:
        chrom, start, end = line.split()[1:4]
        if chrom.replace('chr', '') == args.chromosome_name:
            regions_to_avoid.append([int(start), int(end)])

    regions_to_avoid = sorted(regions_to_avoid, key=itemgetter(0), reverse=False)
    avoid_start = regions_to_avoid[0][0]
    avoid_end   = regions_to_avoid[0][1]
    regions_to_avoid = regions_to_avoid[1:]

    chrom_lens = {}
    with args.chrom_lens_file as chrom_lens_file_lines:
        for line in chrom_lens_file_lines:
            try:
                chrom, chrom_len = line.split()
                chrom_lens[chrom] = int(chrom_len)
            except:
                pass


    chrom_len = chrom_lens[args.chromosome_name]

    start = avoid_end + args.distance_between_SVs
    while SV_lens != []:
        SV = SV_lens[randint(0, len(SV_lens)-1)]
        sv_type, sv_size = SV.split('_')
        end = start + int(sv_size)
        if end > chrom_len:
            print 'Impossible to finish process of creating SVs. Not enough chromosome length (size of %s is %d) for simulated SVs.' % (args.chromosome_name, chrom_len)
            print 'Exiting.'
            break
        
        if between(start, end, avoid_start-args.distance_between_SVs) or between(avoid_start-args.distance_between_SVs, avoid_end+args.distance_between_SVs, start):
            start = avoid_end + int(1.1*args.distance_between_SVs)
            avoid_start = regions_to_avoid[0][0]
            avoid_end   = regions_to_avoid[0][1]
            regions_to_avoid = regions_to_avoid[1:]
        else:
            args.output_bed_file.write('%s\t%d\t%d\t%s\n' % (args.chromosome_name, start, end, sv_type))
            if args.verbose:
                print 'Created ', args.chromosome_name, str(start), str(end), sv_type
            SV_lens.remove(SV)
            start = end + args.distance_between_SVs


if __name__ == '__main__':
    main()
