#include "htslib.hpp"

vcf_parser::vcf_parser(const std::string& filename){
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

int
vcf_parser::open(std::string& filename){
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

int 
vcf_parser::read_snp(){
	int ret = bcf_read(inf, hdr, rec);
	bcf_unpack(rec, BCF_UN_ALL);
	return(ret);
}
bool 
vcf_parser::is_snp(){
	return(bcf_is_snp(rec));
}

bool 
vcf_parser::has_field(std::string field){
	int i,tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, field.c_str());
	return(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id));
}
int 
vcf_parser::extract_field(std::string field){
	nf = bcf_get_format_int32(hdr, rec, field.c_str(), &f, &nf_arr);
	num_entries_per_sample = nf/n;
	return(nf);
}
int 
vcf_parser::extract_GT(){
	ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
	ploidy = ngt/n;
	return(ngt);
}

