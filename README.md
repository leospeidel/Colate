# Colate

**Software for inferring coalescence rates for unphased, low-coverage genomes.**
- Input fileformats: bcf/bam files 
- Requires mutation ages obtained from a genealogy (e.g., using [Relate](https://myersgroup.github.io/relate/)). Precomputed mutation ages available below.

Please send any comments, questions and bug reports to leo.speidel@outlook.com.


# Installation

Precompiled binaries can be found under ./binaries.

Requirements:

- C++11
- cmake
- ZLIB
- OpenSSL
- lzma

To compile from source use
```` bash
cd build
cmake ..
make
````

Colate uses a modifed version of [htslib](https://github.com/samtools/htslib), which has been copied to ./include/vcf/.

# Documentation

## Files required:
- Two sequences, in bcf or bam format, for which coalescence rates will be calculated
- A reference genome in fasta format.
- Genome mask file (optional), to filter out unreliable regions.
- A genealogy to date mutations (see below for precomputed genealogies)

We provide a preprocessed file for the SGDP data here: [SGDP_mutages.tar](https://www.dropbox.com/s/65qbk4lzg50ob34/SGDP_mutages.tar?dl=0) (654Mb).<br/>
Link to human ancestral genomes can be found [here](https://myersgroup.github.io/relate/input_data.html#Data).

## Step 1

**(not needed if using SGDP genealogy from above, you can simply download these files)**

Get mutations ages from a genealogy; this step will add fixed mutations to a *.mut file (see [Relate documentation](https://myersgroup.github.io/relate/getting_started.html#Output) for file format), which is needed for Colate.

- ancestral genome in fasta format
- reference genome in fasta format
- anc/mut files of a Relate genealogy
- bcf of samples used to build Relate genealogy

```` bash
chr=1
${PATH_TO_BINARY}/Colate --mode preprocess_mut \
	--anc example_chr${chr}.anc.gz \
	--mut example_chr${chr}.mut.gz \
	--reference_vcf example_chr${chr}.bcf \
	--ref_genome GRCh37_chr${chr}.fa.gz \
	--anc_genome human_ancestor_${chr}.fa.gz \
	--mask example_mask_chr${chr}.fa \
	-o example_fixed_chr${chr}
````


## Step 2: Run Colate

You can either directly run Colate on bcfs or bams, or precompute an input file, which is beneficial when running Colate repeatedly on the same samples.

### Precompute Colate input files and then run Colate
#### 1. Precompute Colate input files for target and reference: 

**bcfs**:
- bcfs contain samples of interest (Please use e.g., bcftools view -S or -s to subset bcf files).
- chr.txt: Chromosome names, one per line (can be any strings)
- If chromosome names are "1", "2", etc, then input files are *\_chr1.bcf, *\_chr2.bcf etc.
- ref_genome should be separated by chromosome, i.e. GRCh37\_chr1.fa.gz, GRCh37\_chr2.fa.gz etc.
- If --chr is not specified, full file names are required (e.g., --target example.bcf), however these should only contain a single chromosome.

```` bash
mut="example_fixed" #name of .mut files obtained from step 1 (or downloaded)
${PATH_TO_BINARY}/Colate \
	--mode make_tmp \
	--mut ${mut} \
	--target_vcf example_target \ 
	--ref_genome GRCh37 \
        --chr chr.txt \
	-o example_out
````

**bam**:
- bam should be sorted by chromosome name (same order as in chr.txt), and by position.
- ref_genome should be separated by chromosome, i.e. GRCh37\_chr1.fa.gz, GRCh37\_chr2.fa.gz etc.
- chr.txt: Chromosome names, one per line (can be any strings but consistent with file names)

```` bash
mut="example_fixed" #name of .mut files obtained from step 1 (or downloaded)
${PATH_TO_BINARY}/Colate \
	--mode make_tmp \
	--mut ${mut} \
	--target_bam example_target.bam \ 
	--ref_genome GRCh37 \
        --chr chr.txt \
	-o example_out
````

#### 2. Compute coalescenece rates

```` bash
#--target_mask, --reference_mask, 
#--target_age, --reference_age, 
#--years_per_gen, --num_bootstrap are optional

mut="example_fixed" #name of .mut files obtained from step 1 (or downloaded)
bins="3,7,0.2"
${PATH_TO_BINARY}/Colate \
	--mode mut \
	--mut ${mut} \
	--target_tmp example_target.colate.in \
	--reference_tmp example_reference.colate.in \
	--target_mask target_mask \
	--reference_mask reference_mask \
        --bins ${bins} \
	--chr chr.txt \
	--num_bootstrap 20 \
	--target_age 1000 \
	--reference_age 0 \
	--years_per_gen 28 \
	-o example_out
````

### Run directly on bcfs or bams

#### bcfs

- If chromosome names are "1", "2", etc, then input files are *\_chr1.bcf, *\_chr2.bcf etc.
- If --chr is not specified, full file names are required (e.g., --target example.bcf), however these should only contain a single chromosome.
- chr.txt: Chromosome names, one per line (can be any strings)

```` bash
#--target_mask, --reference_mask, 
#--target_age, --reference_age, 
#--years_per_gen, --num_bootstrap are optional
#age is specified in years.
#Assumption is any site not in the bcf is homozygous reference (unless masked out by a mask file).

mut="example_fixed" #name of .mut files obtained from step 1 (or downloaded)
bins="3,7,0.2"
${PATH_TO_BINARY}/Colate \
	--mode mut \
	--mut ${mut} \
	--target_vcf example_target \
	--reference_vcf example_reference \
        --ref_genome GRCh37_chr${chr}.fa.gz \
	--target_mask target_mask \
	--reference_mask reference_mask \
	--bins ${bins} \
	--chr chr.txt \
	--num_bootstrap 20 \
	--target_age 1000 \
	--reference_age 0 \
	--years_per_gen 28 \
	-o example_out
````

#### bam

- bams should be sorted by chromosome name (same order as in chr.txt), and by position.
- ref_genome should be separated by chromosome, i.e. GRCh37\_chr1.fa.gz, GRCh37\_chr2.fa.gz etc.
- chr.txt: Chromosome names, one per line (can be any strings)

```` bash
#--target_age, --reference_age, 
#--years_per_gen, --num_bootstrap are optional

mut="example_fixed" #name of .mut files obtained from step 1 (or downloaded)
bins="3,7,0.2"
${PATH_TO_BINARY}/Colate \
	--mode mut \
	--mut ${mut} \
	--target_bam example_target.bam \
	--reference_bam example_reference.bam \
        --ref_genome GRCh37 \
	--bins ${bins} \
	--chr chr.txt \
	--num_bootstrap 20 \
	--target_age 1000 \
	--reference_age 0 \
	--years_per_gen 28 \
	-o example_out
````


# Output format

*coal files (see [Relate documentation](https://myersgroup.github.io/relate/modules.html#PopulationSizeScript_FileFormats) for file format).
These can be loaded using the [relater R package](https://github.com/leospeidel/relater).

```` R
library(relater)
coal <- read.coal("example_out.coal") #group2 indexes bootstrap iterations
````
Summarize bootstrap iterations:
```` R
library(dplyr)
coal %>% group_by(epoch.start) %>% 
         summarize(mean = mean(haploid.coalescence.rate), 
	           lower = quantile(haploid.coalescence.rate, prob = 0.025), 
		   upper = quantile(haploid.coalescence.rate, prob = 0.975)) -> coal
````

