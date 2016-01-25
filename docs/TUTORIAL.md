# SVGen

LA Lima, H Yang, K Wang. SVGen: Simulation of structural variants in next-generation sequencing data.

<!--
![SVGen workflow](images/svgen_workflow.png?raw=true "SVGen workflow")
-->

## Table of Contents
1. [Download program](#download-program)
2. [Download databases](#databases)
3. [Insert single-nucleotide variants](#insert-SNVs)
4. [Generate structural variant regions](#generate-SVs)
5. [Insert structural variants](#insert-SVs)
6. [Create reads](#create-reads)

## Download program

Using git:

	git clone https://github.com/WGLab/SVGen/
	cd SVGen
	
If you don't use git:

	wget https://github.com/WGLab/SVGen/archive/master.zip
	unzip master.zip
	mv SVGen-master SVGen
	cd SVGen

## Download databases

This step takes some time, but just need to be run once. This command will download the chosen reference, with all chromosomes (from UCSC website) and allele frequencies for different populations (from Annovar website). 

<!-- After downloading, the files will be prepared (indexed...) -->
	

	
	./download_and_format_database.sh $genome_version
	
Example:
	
	
	./download_and_format_database.sh hg38
	
## Insert single-nucleotide variants (SNVs)

Command syntax:	

	python insert_SNVs.py \
        --fasta_input [reference.fa] \
        --fasta_output [output.fa] \
        --freq_file [frequency_file.txt] \
        --chrom [chromosome_name] \
        --vcf_output [output.vcf]
	

Example:

	python insert_SNVs.py \
        --fasta_input reference/hg38/chr22.fa \
        --fasta_output chr22.SNV.fa \
        --freq_file reference/hg38_AFR.sites.2015_08.chrom22.txt \
        --chrom 22 \
        --vcf_output chr22.SNV.vcf


## Generate structural variant (SVs) regions

Command syntax:

    python simulate_SV_BED.py \
        --dup_lens [duplication_lengths_file.txt] \
        --del_lens [deletion_lengths_file.txt] \
        --inv_lens [inversion_lengths_file.txt] \
        --bal_trans_lens [balanced_translocations_lengths_file.txt] \
        --unb_trans_lens [unbalanced_translocations_lengths_file.txt] \
        --chroms [chromosomes_range] \
        --chroms_trans [translocated_chromosomes_range] \
        --chrom_lens [chromosome_lengths.txt] \
        --gaps [genome_gaps.txt] \
        -o [output_file.bed]


Example:

    python simulate_SV_BED.py \
        --dup_lens SV_lengths.txt \
        --del_lens SV_lengths.txt \
        --inv_lens SV_lengths.txt \
        --bal_trans_lens SV_lengths.txt \
        --unb_trans_lens SV_lengths.txt \
        --chroms 22 \
        --chroms_trans 1-10 \
        --chrom_lens reference/chrom_lengths_hg38.txt \
        --gaps reference/gaps_hg38.txt \
        -o SVs.bed
        
## Insert structural variants (SVs)

Command syntax:

    python insert_SVs.py \
        -i [input.fa] \
        -o [output.fa] \
        --chrom_lens [chromosome_lengths.txt] \
        --chrom [chromosome_name] \
        --bed [input_file.bed] \
        -v
        
Example:

    python insert_SVs.py \
        -i chr22.SNV.fa \
        -o chr22.SNV.SV.fa \
        --chrom_lens reference/chrom_lengths_hg38.txt \
        --chrom 22 \
        --bed SVs.bed \
        -v

## Create reads

Command syntax (for **p**aired-**e**nd reads):

    python create_reads.py \
        -pe \
        -i [input.fa] \
        -o [output.fq] \
        --cov [coverage] \
        --read_len [read_length] \
        --snp_rate [SNP_rate_for_reads] \
        --del_rate [deletion_rate_for_reads] \
        --ins_rate [insertion_rate_for_reads] \
        --read_label [label_for_reads] \
        -v
        
Example:

    python create_reads.py \
        -pe \
        -i chr22.SNV.SV.fa \
        -o reads.fq \
        --cov 10 \
        --read_len 100 \
        --snp_rate 0.01 \
        --del_rate 0.0001 \
        --ins_rate 0.0001


	