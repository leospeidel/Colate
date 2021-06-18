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

    int mapq_th = 20;
		int len_th  = 30;
		int bq_th   = 30;
		int mismatch_th = 3;

    void count_alleles_for_read();

  public:

    samFile *fp_in; //bam file
    bam_hdr_t *bamHdr; //header
    bam1_t *aln; //alignment

    bool eof;

    int32_t pos; //left most position of alignment in zero-based coordinate
    char *chr = NULL; //contig name
    uint32_t len; //length of the read
    uint32_t maxlen;
    uint8_t *q; //quality string
    uint32_t mapq; //mapping quality
    char* seq; //seq

    std::string contig = "";
    fasta ref_genome;

    double coverage = 0, coverage_after_filter = 0;

    int num_entries = 1e5; //num_entries in pos_of_entry/count_alleles
    std::vector<int> pos_of_entry;
    std::vector<std::vector<int>> count_alleles;

    bam_parser();
    bam_parser(const std::string& filename, const std::string& params);
    bam_parser(const std::string& filename, const std::string& params, const std::string& filename_ref);
    ~bam_parser(){
      bam_destroy1(aln);
      sam_close(fp_in);
    }

    int open(std::string& filename, std::string& filename_ref);
    int read_entry();
    //read all reads covering current_pos
    bool read_to_pos(int current_pos);
    void assign_contig(std::string& icontig, std::string& filename_ref);

};

#endif //HTSLIB_HPP

