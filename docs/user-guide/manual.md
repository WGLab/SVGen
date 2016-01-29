# Commands parameters


## Table of Contents
1. [Download databases](#databases)
2. [Insert single-nucleotide variants](#insert-snv)
3. [Simulate structural variant regions](#simulate-sv)
4. [Insert structural variants](#insert-sv)
5. [Create reads](#create-reads)


### Download databases
	
Command syntax:
	
	./download_and_format_database.sh $genome_version
	
There are only two ways to run this command:
	
	./download_and_format_database.sh hg19
	
or
		
	./download_and_format_database.sh hg38
	
## Insert single-nucleotide variants
This command inserts single-nucleotide variants (SNVs) in a given fasta file.

	python insert_SNVs.py \
        --fasta_input [reference.fa] \
        --fasta_output [output.fa] \
        --freq_file [frequency_file.txt] \
        --chrom [chromosome_name] \
        --vcf_output [output.vcf]


Required arguments:
	  
	  --fasta_input [input.fasta]
                           Input fasta file to be used as reference to
                           receive single-nucleotide variants (SNVs).
	                        
	  --fasta_output [output.fasta]
                           Output fasta file with random SNVs,
                           based in frequencies.
	                        
	  --freq_file [freq.txt]
                           Text file SNV frequencies.
	                        
	  --vcf_output [output.vcf]
                           VCF file generated with SNVs inserted.
	                        
	  --chrom [chromosome]
                           Chromosome name.


Optional arguments:

	  -h/--help            
                           Shows the help message and exits


## Simulate structural variant regions

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

Required arguments:

	  --chrom_lens [chrom_lens.txt]
                           Text file with chromosome lengths.
	                        
	  --output/-o [output.bed]
                           BED output file.

	  --gaps [gaps_file.bed]      
                           BED file with regions to avoid
                           (centromeres and telomeres).
                           
	  --chroms [chrom_range]
                           Chromosome names to generate SV regions.
                           (e.g. 3, 1-5, 20-Y)
                           (default: all)

Optional arguments:

	  -h, --help            
	                        Shows the help message and exits
	                        
	  --del_lens [del_lens.txt]
                           Text file with deletion
                           lengths (one per line).
	                        
	  --dup_lens [dup_lens.txt]
                           Text file with duplication
                           lengths (one per line).
	                        
	  --inv_lens [inv_lens.txt]
                           Text file with inversion
                           lengths (one per line).
	                        
	  --bal_trans_lens [bal_lens.txt]
                           Text file with balanced translocation
                           lengths (one per line).
	                        
	  --unb_trans_lens [unb_lens.txt]
                           Text file with unbalanced translocation
                           lengths (one per line). 
	                        
	  --chroms_trans [chrom_range]
                           Chromosomes from which translocations
                           will come from (e.g. 3, 1-5, 20-Y)
                           (default: all)
	                        
	  --distance/-d [distance between SVs]
                           Average distance between SVs in a
                           countinuous (ungapped) region.
                           (default: 100000)
                           
	  --dist_sd/-sd [distance st. dev.]
                           Standard deviation of distance between
                           SVs in a countinuous (ungapped) region.
                           (default: 10000)
                                                      
	  -v/--verbose
                           Shows more details about program execution.
     
## Insert structural variants


    python insert_SVs.py \
        -i [input.fa] \
        -o [output.fa] \
        --chrom_lens [chromosome_lengths.txt] \
        --chrom [chromosome_name] \
        --bed [SVs.bed] \
        -v

Required arguments:

	  
     --fasta_input [input.fasta]
                           Input fasta file to be used as reference to
                           receive structural variants (SVs).
	                        
     --fasta_output [output.fasta]
                           Output fasta file with SVs,
                           based on BED file.

     --chrom_lens [chrom_lens.txt]
                           Text file with chromosome lengths.

     --bed [SVs.bed]
                           BED file with SVs to be inserted
                           (it can be provided by user, or just
                           the output of 'simulate_SV_BED.py').

Optional arguments:

     --fasta_label [fasta_label]
                           Name to label fasta sequences.

     -h/--help            
                           Shows the help message and exits.
	                        
     -v/--verbose
                           Shows more details about program execution.
	                        

## Create reads

Example with paired-end reads:

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
        

Required arguments:
	
     --fasta_input/-i [input.fasta]
                           Fasta file to be changed with SVs.

     --output_file/-o [output.fq]
                           Output file for reads. It must finish with
                           .fq/.fastq (paired-end option will automatically
                           change file names to output1.fq and output2.fq).

     --cov [coverage]      
                           Average coverage.
	  
     --read_len [avg read len]
                           Average read length (this length will be fix
                           for short reads, but not for long reads).
                           This parameters will determine whether to
                           create reads according to Illumina error rate
                           (length <= 500) or PacBio error rate (length > 500).
	                        

	  
Optional arguments:

     -h/--help            
                           Shows the help message and exits.

     -pe                   
                           Add option to generate paired-end reads.
                           If not given, single-end reads will be
                           generated.
                           
     --ins_rate [insertion rate]
                           Insertion error rate for reads.
                           (default: 0.12)
	                        
     --del_rate [deletion rate]
                           Deletion error rate for reads.
                           (default: 0.02)	                        
     --snp_rate [snp rate]
                           SNP error rate for reads.
                           (default: 0.01)
	  
     --insert_size [insert size]
                           Insert size for short reads.
                           (default: 300)	                        
     --insert_sd [insert st. dev.]
                           Insert standard deviation for short reads.
                           (default: 50)
	                        
     --alpha [alpha value]
                           Alpha value for beta distribution of base qualities.
                           (default: 2)
	  
     --beta [beta value]           
                           Beta for beta distribution of read lengths.
                           (default: 10)
	  
	  --read_label [label]
	                        Label to add in each read.
	                        
      -v/--verbose
                           Shows more details about program execution.
	                        

	
