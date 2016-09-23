# SVGen: Simulation of structural variants in next-generation sequencing data.

Existing simulation tools for structural variants (SV) offer limited options, some do not simulate single-nucleotide variants (SNV), and all need external programs to generate simulated sequence reads to benchmark SV calling software. To address these limitations, we developed SVGen to be a tool with an extensive range of options to simulate germline and somatic SVs of various types and from short and long read sequencing platforms. The output from SVGen include BED files with SV regions, FASTA files with the simulated genome sequence and FASTQ files for short or long reads.

## Features

* Insertion of single-nucleotide variations (SNVs), according to a frequencies file (provided for different populations).

* Simulation of different types of structural variations (SVs): deletions, duplications, inversions, balanced and unbalanced translocations.

* Simulation of short (single- and paired-end) and long reads, following error models from Illumina and PacBio reads, respectively.

* Easy download and no need of installation. No dependencies/requirements, just Python 2.7.


## Synopsis

* **download and prepare databases**: download databases and allele frequencies files from UCSC and Annovar website

* **insert SNVs**: from an allele frequency file, insert SNVs in a given FASTA file.

* **simulated SV regions**: creates a BED file with up to 5 types of structural variant regions: deletions, duplications, inversions, balanced and unbalanced translocations.

* **insert SVs**: given a BED file with SVs and a FASTA file, inserts the correponding SVs into the FASTA file.

* **simulated reads**: generates short or long reads from a given FASTA file.

## Tutorial

Click [here](http://svgen.openbioinformatics.org/en/latest/user-guide/startup/) for a **quick start**, or [here](http://svgen.openbioinformatics.org/en/latest/user-guide/manual/) for a **detailed tutorial**.


## Revision History

For details, please go [here](https://github.com/WGLab/SVGen/commits/master).

## Contact

Leandro Lima - leandrol@usc.edu / lelimaufc@gmail.com



<!-- Please join [SVGen users group](https://groups.google.com/forum/#!forum/svgen-users) for updates! -->

## Citation

LA Lima, H Yang, C Dong, K Wang. SVGen: Simulation of structural variants in next-generation sequencing data. (submitted)

## License

[WGLab MIT License](http://wglab.mit-license.org)

## More information

* [SVGen Homepage](http://svgen.openbioinformatics.org)

* [Wang Genomics Lab Homepage](http://wglab.org)



