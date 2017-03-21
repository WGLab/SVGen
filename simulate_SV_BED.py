#!/usr/bin/python2
# simulate_SV_BED.py
# C: Oct  1, 2015
# M: Jan 28, 2016
# A: Leandro Lima <leandrol@usc.edu>


import sys, argparse
from random import randint, choice, gauss
from operator import itemgetter

prog_name = 'simulate_SV_BED.py'

all_chroms = map(str, range(1,23)) + ['X', 'Y']

def pick_an_interval(size, intervals, distance):
    interval = choice(intervals)
    intervals.remove(interval)
    # Looking for an interval that has enough room for the chosen size
    while interval[1] - interval[0] + 1 < distance + size and intervals != []:
        interval = choice(intervals)
        intervals.remove(interval)
    if interval[1] - interval[0] + 1 < size + distance:
        print '\n\tError! It was not possible to finish creating regions. Try to run again with less or smaller regions, or a smaller distance.\n'
        return None, None, intervals
    start = interval[0] + distance
    end   = start + size - 1
    interval[0] = end + 1
    if interval[0] < interval[1]:
        intervals.append(interval)

    return start, end, intervals


def intervals_overlap(start_1, end_1, start_2, end_2):
    min_1 = min(start_1, end_1)
    max_1 = max(start_1, end_1)
    min_2 = min(start_2, end_2)
    max_2 = max(start_2, end_2)
    if max_2 < min_1 or max_1 < min_2: # no overlap
        return False
    else:
        return True


"""
def between(start, end, value):
    if (start <= value and value <= end) or (end <= value and value <= start):
        return True
    else:
        return False
"""


def remove_intervals(intervals, start, end):
    new_intervals = []
    intervals_sorted = sorted(intervals, key=lambda tup: tup[0])
    interval_start = min(start, intervals_sorted[0][1] + 1)
    for i in range(len(intervals_sorted)):
        interval_end   = intervals_sorted[i][0] - 1
        if i > 0 or interval_start < interval_end:
            new_intervals.append([interval_start, interval_end])
        interval_start = intervals_sorted[i][1] + 1
    if intervals_sorted[-1][1]+1 < end:
        new_intervals.append([intervals_sorted[-1][1]+1, end])
    return new_intervals



# Function from https://codereview.stackexchange.com/a/69249/89845
def sort_and_merge_intervals(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)
    return merged


"""
def pick_a_region(size, start, regions_to_avoid, distance):
    avoid_start = regions_to_avoid[0][0]
    avoid_end   = regions_to_avoid[0][1]

    end = start + size - 1

    while avoid_start < start:
        regions_to_avoid = regions_to_avoid[1:]
        avoid_start = regions_to_avoid[0][0]
        avoid_end   = regions_to_avoid[0][1]


    # while between(avoid_start, avoid_end, start) or between(avoid_start, avoid_end, end):
    while intervals_overlap(avoid_start, avoid_end, start, end):
        print 'jump ',
        start = avoid_end + distance
        end = start + size - 1
        regions_to_avoid = regions_to_avoid[1:]
        avoid_start = regions_to_avoid[0][0]
        avoid_end   = regions_to_avoid[0][1]
       
    print 
    print "pick_a_region"
    print avoid_start, avoid_end, start, end, regions_to_avoid[:3]
    return start, end, regions_to_avoid
"""


def parse_chrom_ranges(text):
    chroms = []
    chroms_splitted = text.replace('chr', '').split(',')
    if chroms_splitted == ['all']:
        chroms = all_chroms
    else:
        for piece in chroms_splitted:
            piece = piece.replace('X', '23')
            piece = piece.replace('Y', '24')
            try:
                first, last = piece.split('-')
                for c in range(int(first), int(last)+1):
                    chroms.append(str(c))
            except:
                chroms.append(piece)
    if '23' in chroms:
        chroms.remove('23')
        chroms.append('X')
    if '24' in chroms:
        chroms.remove('24')
        chroms.append('Y')
    return list(set(chroms))



def main():

    parser = argparse.ArgumentParser(description='Get arguments to create random structural variant regions in a BED file.', prog=prog_name)
    # Required
    parser.add_argument('--chrom_lens', required=True, type=file, metavar='chrom_lengths_file', dest='chrom_lens_file', help='Text file with chromosome lengths.')
    parser.add_argument('--output', '-o', required=True, type=argparse.FileType('w'), metavar='output_bed_file', dest='output_bed_file', help='BED output file.')
    parser.add_argument('--gaps', required=True, type=file, metavar='gaps_file', dest='avoid_regions_file', help='BED file with regions to avoid (centromeres and telomeres).')
    parser.add_argument('--chroms', required=True, type=str, metavar='chromosome_names', dest='chromosome_name', help='Chromosome names (range).', default='all')
    # Optional
    parser.add_argument('--del_lens', required=False, type=file, metavar='del_lengths_file', dest='del_len_filename',   help='Text file with deletion lengths.')
    parser.add_argument('--dup_lens', required=False, type=file, metavar='dup_lengths_file', dest='dup_len_filename',   help='Text file with duplication lengths.')
    parser.add_argument('--inv_lens', required=False, type=file, metavar='inv_lengths_file', dest='inv_len_filename',   help='Text file with inversion lengths.')
    parser.add_argument('--bal_trans_lens', required=False, type=file, metavar='bal_trans_lengths_file', dest='bal_trans_len_filename', help='Text file with balanced translocation lengths.')
    parser.add_argument('--unb_trans_lens', required=False, type=file, metavar='unb_trans_lengths_file', dest='unb_trans_len_filename', help='Text file with unbalanced translocation lengths.')
    parser.add_argument('--chroms_trans', required=False, type=str, metavar='chroms_trans', help='Chromosomes from which translocations will come from.', default='all')
    parser.add_argument('--distance', '-d', required=False, type=int, metavar='distance_between_SVs', dest='distance_between_SVs', help='Distance between SVs in a countinuous (ungapped) region.', default=100000)
    parser.add_argument('--dist_sd', '-sd', required=False, type=int, metavar='distance_sd', dest='distance_sd', help='Standard deviation of distance between SVs in a countinuous (ungapped) region.', default=10000)

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()


    # Checking if at least one type of variant was provided
    if args.del_len_filename is None and \
       args.dup_len_filename is None and \
       args.inv_len_filename is None and \
       args.bal_trans_len_filename is None and \
       args.unb_trans_len_filename is None:
        print '\n\n\tError! At least one file with variant (deletion, duplication, inversion or translocation) lengths has to be provided.\n\n'
        sys.exit(1)


    chrom_lens = {}
    with args.chrom_lens_file as chrom_lens_file_lines:
        for line in chrom_lens_file_lines:
            try:
                chrom, chrom_len = line.split()
                chrom_lens[chrom] = int(chrom_len)
            except:
                pass


    # Avoiding gap regions (centromeres, telomeres, etc.): https://www.biostars.org/p/2349/#2351
    lines = args.avoid_regions_file.read().split('\n')
    while lines[-1] == '':
        lines.pop()

    while lines[0].startswith('#'):
        lines = lines[1:]

    regions_to_avoid = {}
    nongap_regions = {}
    for line in lines:
        chrom, start, end = line.split()[1:4]
        chrom = chrom.replace('chr', '')
        if not chrom in regions_to_avoid.keys():
            regions_to_avoid[chrom] = [[int(start), int(end)]]
        else:
            regions_to_avoid[chrom] += [[int(start), int(end)]]


    # Parse chromosome ranges to get individual chromosomes
    chroms = parse_chrom_ranges(args.chromosome_name)
    chroms_trans = parse_chrom_ranges(args.chroms_trans)


    for chrom in all_chroms:
        gap_regions = sort_and_merge_intervals(regions_to_avoid[chrom])
        nongap_regions[chrom] = remove_intervals(gap_regions, 1, chrom_lens[chrom])
        # print
        # print 'chrom', chrom, '>>', chrom_lens[chrom]
        # print 'gap_regions', gap_regions
        # print 'nongap_regions', nongap_regions[chrom]


    SV_codes = {}
    SV_regions = {}
    del_lens = []
    dup_lens = []
    inv_lens = []
    bal_trans_lens = []
    unb_trans_lens = []

    for chrom in parse_chrom_ranges('all'):
        SV_regions[chrom] = []

    
    for chromA in chroms:
        if args.verbose:
            print 'Creating SVs for chromosome', chromA
        SV_codes[chromA] = []

        if not args.del_len_filename is None:
            if del_lens == []:
                lines = args.del_len_filename.read().split('\n')
                while lines[-1] == '':
                    lines.pop()
                del_lens = lines
            SV_codes[chromA] += ['del_'+length for length in del_lens]

        if not args.dup_len_filename is None:
            if dup_lens == []:
                lines = args.dup_len_filename.read().split('\n')
                while lines[-1] == '':
                    lines.pop()
                dup_lens = lines
            SV_codes[chromA] += ['dup_'+length for length in dup_lens]

        if not args.inv_len_filename is None:
            if inv_lens == []:
                lines = args.inv_len_filename.read().split('\n')
                while lines[-1] == '':
                    lines.pop()
                inv_lens = lines
            SV_codes[chromA] += ['inv_'+length for length in inv_lens]

        if not args.bal_trans_len_filename is None:
            if bal_trans_lens == []:
                lines = args.bal_trans_len_filename.read().split('\n')
                while lines[-1] == '':
                    lines.pop()
                bal_trans_lens = lines
            for length in bal_trans_lens:
                chromB = choice(chroms_trans)
                SV_codes[chromA] += ['baltr_'+length+'_'+chromB]

        if not args.unb_trans_len_filename is None:
            if unb_trans_lens == []:
                lines = args.unb_trans_len_filename.read().split('\n')
                while lines[-1] == '':
                    lines.pop()
                unb_trans_lens = lines
            for length in unb_trans_lens:
                chromB = choice(chroms_trans)
                SV_codes[chromA] += ['unbtr_'+length+'_'+chromB]


        chromA_len = chrom_lens[chromA]

        endA = 1

        while SV_codes[chromA] != []:
            chromB = None
            SV = choice(SV_codes[chromA])
            try:
                sv_type, sv_size = SV.split('_')
            except:
                sv_type, sv_size, chromB = SV.split('_')

            sv_size = int(sv_size)

            # startA, endA, regions_to_avoid[chromA] = pick_a_region(sv_size, endA + int(gauss(args.distance_between_SVs, args.distance_sd)), regions_to_avoid[chromA], int(gauss(args.distance_between_SVs, args.distance_sd)))
            startA, endA, nongap_regions[chromA] = pick_an_interval(sv_size, nongap_regions[chromA], int(gauss(args.distance_between_SVs, args.distance_sd)))

            if endA > chromA_len:
                print '\n\n\tImpossible to finish process of creating SVs for chrom %s. Not enough chromosome length (size is %d) for simulated SVs.' % (chromA, chromA_len)
                print '\tExiting.\n\n'
                break
            
            if chromB is None:
                SV_regions[chromA] += [sv_type + '_' + str(startA) + '_' + str(endA)]
            else:
                startB, endB, nongap_regions[chromB] = pick_an_interval(sv_size, nongap_regions[chromB], int(gauss(args.distance_between_SVs, args.distance_sd)))
                '''
                if SV_regions[chromB] == []:
                    # If no SV was inserted in chromB, start from the first position in chromosome
                    # startB, endB, regions_to_avoid[chromB] = pick_a_region(sv_size, 1, regions_to_avoid[chromB], int(gauss(args.distance_between_SVs, args.distance_sd)))
                    startB, endB, nongap_regions[chromB] = pick_an_interval(sv_size, nongap_regions[chromA], int(gauss(args.distance_between_SVs, args.distance_sd)))
                else:
                    # Else, get the last SV inserted and start from the end of such SV
                    sv_typeB, startB, endB = SV_regions[chromB][-1].split('_')[:3]
                    startB, endB, regions_to_avoid[chromB] = pick_a_region(sv_size, int(endB) + int(gauss(args.distance_between_SVs, args.distance_sd)), regions_to_avoid[chromB], int(gauss(args.distance_between_SVs, args.distance_sd)))
                    # If chromosomes are the same, pick another region
                    if chromA == chromB:
                        startB, endB, regions_to_avoid[chromB] = pick_a_region(sv_size, int(endB) + int(gauss(args.distance_between_SVs, args.distance_sd)), regions_to_avoid[chromB], int(gauss(args.distance_between_SVs, args.distance_sd)))
                '''

                SV_regions[chromA] += [sv_type + '_' + str(startA) + '_' + str(endA) + '_' + chromB + ':' + str(startB) + '-' + str(endB)]
                if sv_type == 'unbtr':
                    SV_regions[chromB] += ['unaff_' + str(startB) + '_' + str(endB) + '_' + chromA + ':' + str(startA) + '-' + str(endA)]
                else:
                    SV_regions[chromB] += [sv_type + '_' + str(startB) + '_' + str(endB) + '_' + chromA + ':' + str(startA) + '-' + str(endA)]


            SV_codes[chromA].remove(SV)

    for chrom in all_chroms:
        for SV in SV_regions[chrom]:
            elems = SV.split('_')
            sv_type = elems[0]
            if sv_type != 'unaff': # This means that the region in the current chromosome was not affected
                start   = elems[1]
                end     = elems[2]
                if len(elems) > 3:
                    info = elems[3]
                    args.output_bed_file.write('%s\t%s\t%s\t%s\t%s\n' % (chrom, start, end, sv_type, info))
                else:
                    args.output_bed_file.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, sv_type))


if __name__ == '__main__':
    main()
