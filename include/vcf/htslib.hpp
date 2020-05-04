#include <stdio.h>
#include <iostream>
#include <string>


#include "htslib/vcf.h"

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
		vcf_parser(const std::string& filename){
			inf = bcf_open(filename.c_str(), "r");
			if (inf == NULL) {
			  std::cerr << "Warning: Failed to open file " << filename << std::endl;
			}else{
				hdr = bcf_hdr_read(inf);
				// struc for storing each record
				rec = bcf_init();
				n = bcf_hdr_nsamples(hdr);
			}
		}
		~vcf_parser(){
			free(f);
			bcf_hdr_destroy(hdr);
			bcf_close(inf);
			bcf_destroy(rec);
		}

		int open(std::string& filename){
			inf = bcf_open(filename.c_str(), "r");
			if (inf == NULL) {
				return EXIT_FAILURE;
			}
			hdr = bcf_hdr_read(inf);
			// struc for storing each record
			rec = bcf_init();
			n = bcf_hdr_nsamples(hdr);
			return EXIT_SUCCESS;
		}

		int read_snp(){
			int ret = bcf_read(inf, hdr, rec);
			bcf_unpack(rec, BCF_UN_ALL);
			return(ret);
		}
		bool is_snp(){
			return(bcf_is_snp(rec));
		}

		bool has_field(std::string field){
			int i,tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, field.c_str());
			return(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id));
		}
		int extract_field(std::string field){
			nf = bcf_get_format_int32(hdr, rec, field.c_str(), &f, &nf_arr);
			num_entries_per_sample = nf/n;
			return(nf);
		}
		int extract_GT(){
			ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
			ploidy = ngt/n;
			return(ngt);
		}


};
