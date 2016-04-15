1. **How does SVGen differ from other simulation tools?**

    SVGen is designed to simulate raw reads (FASTQ files) for benchmarking software tools that detect strcutural variants from sequencing data, yet many other tools aim to simulate regions (BED files). Currently SVGen supports both short reads and long reads (PacBio-generated reads), each with different error models. Furthermore, SVGen can simulate loss-of-heterozygosity (LOH) regions, and simulate SVs in tumor genomes, given a stromal contamination parameter by the user.

