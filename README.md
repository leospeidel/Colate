# Colate

# Requirements

- C++11
- cmake
- ZLIB
- OpenSSL

# Installation

Precompiled binaries can be found under ./binaries.

To compile from source use
```` bash
cd build
cmake ..
make
````

# Documentation

- this will be refined in the coming days.

## Step 1

Get mutations ages from a genealogy.
```` bash
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



## Step 2

Run Colate




## Alternative Step 2

Precompute Colate input files

## Alternative Step 3



