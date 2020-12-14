# Colate

Please send any questions and bug reports to leo.speidel@outlook.com.


# Installation

Precompiled binaries can be found under ./binaries.

Requirements:

- C++11
- cmake
- ZLIB
- OpenSSL

To compile from source use
```` bash
cd build
cmake ..
make
````

# Documentation

- Preliminary - will be refined in future. 

## Step 1

Get mutations ages from a genealogy; this step will add fixed mutations to a *.mut file (see [Relate documentation](https://myersgroup.github.io/relate/getting_started.html#Output) for file format), which is needed for Colate.
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

We provide a preprocessed file for the SGDP data here:

## Step 2: Run Colate

- chr.txt: Chromosome names, one per line (can be any strings)

You can either directly run Colate on vcfs or bams, or precompute an input file, which is beneficial when running Colate repeatedly on the same samples.

### Precompute Colate input files and then run Colate
#### 1. Precompute Colate input files for target and reference: 

**bcfs**:
- If chromosome names are "1", "2", etc, then input files are *\_chr1.bcf, *\_chr2.bcf etc.
- If --chr is not specified, full file names are required (e.g., --target example.bcf), however these should only contain a single chromosome.

```` bash
mut="example_fixed"
${PATH_TO_BINARY}/Colate \
	--mode make_tmp \
	--mut ${mut} \
	--target_vcf example_target \ 
	--ref_genome GRCh37 \
        --chr chr.txt \
	-o example_out
````

**bams**:
- bams should be sorted by chromosome name (same order as in chr.txt), and by position.
- ref_genome should be separated by chromosome, i.e. GRCh37\_chr1.fa.gz, GRCh37\_chr2.fa.gz etc.

```` bash
mut="example_fixed"
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

mut="example_fixed"
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

### Run directory on bcfs or bams 

#### bcfs

- If chromosome names are "1", "2", etc, then input files are *\_chr1.bcf, *\_chr2.bcf etc.
- If --chr is not specified, full file names are required (e.g., --target example.bcf), however these should only contain a single chromosome.

```` bash
#--target_mask, --reference_mask, 
#--target_age, --reference_age, 
#--years_per_gen, --num_bootstrap are optional
#age is specified in years.
#Assumption is any site not in the bcf is homozygous reference (unless masked out by a mask file).

mut="example_fixed"
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

#### bams

- bams should be sorted by chromosome name (same order as in chr.txt), and by position.
- ref_genome should be separated by chromosome, i.e. GRCh37\_chr1.fa.gz, GRCh37\_chr2.fa.gz etc.

```` bash
#--target_age, --reference_age, 
#--years_per_gen, --num_bootstrap are optional

mut="example_fixed"
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
