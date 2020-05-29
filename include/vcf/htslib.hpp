#ifndef HTSLIB_HPP
#define HTSLIB_HPP

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

#include "data.hpp"

#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/sam.h"

class vcf_parser{

	private:

	public:

		int n    = 0;  // total number of records in file

		//stores genotype
		int ngt_arr = 0;
		int ngt     = 0;
		int *gt     = NULL;
		int ploidy = 0;

		//stores any format field
		int nf_arr = 0;
		int nf     = 0;
		int *f     = NULL;
		int num_entries_per_sample = 0;

		htsFile * inf; //file
		bcf_hdr_t * hdr; //header
		bcf1_t * rec; //record

		vcf_parser();
		vcf_parser(const std::string& filename);
		~vcf_parser(){
			free(f);
			bcf_hdr_destroy(hdr);
			bcf_close(inf);
			bcf_destroy(rec);
		}

		int open(std::string& filename);
		int read_snp();
		bool is_snp();
		bool has_field(std::string field);
		int extract_field(std::string field);
		int extract_GT();

};

class bam_parser{

	private:

		int prev_pos;

		int count_alleles_for_read(){

			if(mapq > 30 && len >= 10){

				//calculate how many bp match reference at this read
				int num_matching = 0;
				for(int i = 0; i < len; i++){

					if(ref_genome.seq[pos+i] == seq[i]) num_matching++;

					//reset count_alleles if necessary
					if(pos_of_entry[(pos+i) % num_entries] != pos){
						std::fill(count_alleles[(pos+i) % num_entries].begin(), count_alleles[(pos+i) % num_entries].end(), 0);
						pos_of_entry[(pos+i) % num_entries] = pos+i;
					}

				}
				if(num_matching >= 0.9*len){

					coverage_after_filter += len;
					for(int i = 0; i < len; i++){

						if(seq[i] == 'A'){
							count_alleles[(pos+i) % num_entries][0]++;
						}else if(seq[i] == 'C'){
							count_alleles[(pos+i) % num_entries][1]++;
						}else if(seq[i] == 'G'){
							count_alleles[(pos+i) % num_entries][2]++;
						}else if(seq[i] == 'T'){
							count_alleles[(pos+i) % num_entries][3]++;
						}else{
							assert(1);
						}

					}
				}
			}

		}


	public:

		samFile *fp_in; //bam file
		bam_hdr_t *bamHdr; //header
		bam1_t *aln; //alignment

    bool eof;

		int32_t pos; //left most position of alignment in zero-based coordinate
		char *chr; //contif name
		uint32_t len; //length of the read
		uint8_t *q; //quality string
		uint32_t mapq; //mapping quality
		char* seq; //seq

		fasta ref_genome;

		int coverage = 0, coverage_after_filter = 0;

		int num_entries = 10000; //num_entries in pos_of_entry/count_alleles
		std::vector<int> pos_of_entry;
		std::vector<std::vector<int>> count_alleles;

		bam_parser();
		bam_parser(const std::string& filename, const std::string& filename_ref){
			fp_in = hts_open(filename.c_str(), "r");
			if (fp_in == NULL) {
				std::cerr << "Warning: Failed to open file " << filename << std::endl;
			}else{
				bamHdr = sam_hdr_read(fp_in); //read header
				aln = bam_init1(); //initialize an alignment

				eof = false;
				coverage = 0;
				coverage_after_filter = 0;

				//initialise count_alleles at any given position
				prev_pos = -1;
				pos = 0;
				pos_of_entry.resize(num_entries);
				std::fill(pos_of_entry.begin(), pos_of_entry.end(), -1);
				count_alleles.resize(num_entries);
				std::vector<std::vector<int>>::iterator it_count_alleles;
				for(it_count_alleles = count_alleles.begin(); it_count_alleles != count_alleles.end(); it_count_alleles++){
					(*it_count_alleles).resize(4);
				}

				ref_genome.Read(filename_ref);
			}
		}
		~bam_parser(){
			bam_destroy1(aln);
			sam_close(fp_in);
		}

		int open(std::string& filename, std::string& filename_ref){
			fp_in = hts_open(filename.c_str(), "r");
			if (fp_in == NULL) {
				return EXIT_FAILURE;
			}else{
				bamHdr = sam_hdr_read(fp_in); //read header
				aln = bam_init1(); //initialize an alignment

				eof = false;
				coverage = 0;
				coverage_after_filter = 0;

				//initialise count_alleles at any given position
				prev_pos = -1;
				pos = 0;
				pos_of_entry.resize(num_entries);
				std::fill(pos_of_entry.begin(), pos_of_entry.end(), -1);
				count_alleles.resize(num_entries);
				std::vector<std::vector<int>>::iterator it_count_alleles;
				for(it_count_alleles = count_alleles.begin(); it_count_alleles != count_alleles.end(); it_count_alleles++){
					(*it_count_alleles).resize(4);
				}
				ref_genome.Read(filename_ref);

			}
			return EXIT_SUCCESS;
		}

		int read_entry(){
			int ret = sam_read1(fp_in, bamHdr, aln);

			if(ret > 0){
				pos = aln->core.pos; //left most position of alignment in zero based coordianate (+1)
				chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
				len = aln->core.l_qseq; //length of the read.

				q = bam_get_seq(aln); //quality string
				mapq = aln->core.qual ; //mapping quality

				seq = (char *)malloc(len);
				for(int i=0; i< len ; i++){
					seq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
				}

				count_alleles_for_read();

				if(pos < prev_pos){
					std::cerr << "Error: BAM file not sorted by position." << std::endl;
					exit(1);
				}
				prev_pos = pos;

				coverage += len;
			}else{
        eof = true;
			}
			return(ret);
		}

		//read all reads covering current_pos
		int read_to_pos(int current_pos){
			if( !eof && pos - current_pos < num_entries/2.0){
				while(read_entry() > 0){
					if(pos - current_pos >= num_entries/2.0) break;
				}
			}
		} 

};

#endif //HTSLIB_HPP

