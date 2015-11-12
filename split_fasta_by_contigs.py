# split_fasta_by_contigs.py
# C: Sep 23, 2014
# M: Nov 12, 2015
# A: Leandro Lima <leandrol@usc.edu>

import sys

def main():
    filename = sys.argv[1]
    print 'Opening input file.'    
    fasta = open(filename).read()
    print 'Splitting input file.'
    seqs = fasta.split('>')[1:]
    for seq in seqs:
        label = seq.split('\n')[0].split()[0]
        print 'Creating file %s.fa' % label
        output = open(label + '.fa', 'w')
        output.write('>%s' % seq)
        output.close()
    

if __name__ == '__main__':
    main()
