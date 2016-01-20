# SVGen

LA Lima, H Yang, K Wang. SVGen: Simulation of structural variants in next-generation sequencing data.

<!--
![SVGen workflow](images/svgen_workflow.png?raw=true "SVGen workflow")
-->

## Table of Contents
1. [Download program](#download-program)
2. [Download databases](#databases)
3. [Run](#run)
4. [Other examples](#other-examples)

## Download program

Using git:

	git clone https://github.com/WGLab/SVGen/
	cd SVGen
	
If you don't use git:

	wget https://github.com/WGLab/SVGen/archive/master.zip
	unzip master.zip
	mv SVGen-master SVGen
	cd SVGen

## Download databases (genome references and frequency files).

This step takes some time, but just need to be run once. This command will download the chosen reference, with all chromosomes () and allele frequencies for different populations (from Annovar website). After downloading, the files will be prepared (indexed
	

	
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



	