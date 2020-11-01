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



void 
bam_parser::count_alleles_for_read(){

  if(mapq > 30 && len >= 30){

    //calculate how many bp match reference at this read
    int num_matching = 0, total = 0;
    int start = std::max(5, (int)(len*0.05)), end = len - std::max(5, (int)(len*0.05));

    for(int i = start; i < end; i++){

      total++;
      if(ref_genome.seq[pos+i] == seq[i]) num_matching++;

      //reset count_alleles if necessary
      if(pos_of_entry[(pos+i) % num_entries] != pos+i){
        std::fill(count_alleles[(pos+i) % num_entries].begin(), count_alleles[(pos+i) % num_entries].end(), 0);
        pos_of_entry[(pos+i) % num_entries] = pos+i;
      }

    }

    if(total - num_matching <= 3){

      coverage_after_filter += len;
      for(int i = start; i < end; i++){

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

bam_parser::bam_parser(const std::string& filename){
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
  }
}

bam_parser::bam_parser(const std::string& filename, const std::string& filename_ref){
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
    read_entry();
    contig = chr;
  }
}

int 
bam_parser::open(std::string& filename, std::string& filename_ref){
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

int 
bam_parser::read_entry(){
  int ret = sam_read1(fp_in, bamHdr, aln);

  if(ret > 0){
    pos = aln->core.pos; //left most position of alignment in zero based coordianate (+1)
    chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
    if(strcmp(chr, contig.c_str()) == 0){
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
    }
  }else{
    eof = true;
  }
  return(ret);
}

bool 
bam_parser::read_to_pos(int current_pos){
  if(strcmp(contig.c_str(), chr) == 0){
    if( !eof && pos - current_pos < num_entries/2.0){
      while(read_entry() > 0){
        if(pos - current_pos >= num_entries/2.0) break;
        if(strcmp(contig.c_str(), chr) != 0) break;
      }
    }
  }
  return(eof);
}

void
bam_parser::assign_contig(std::string& icontig, std::string& filename_ref){

  if(icontig != ""){
    contig = icontig;
  }
	ref_genome.seq.clear();
  ref_genome.Read(filename_ref);
  eof = false;
  coverage = 0;
  coverage_after_filter = 0;

  //initialise count_alleles at any given position
  prev_pos = -1;
  pos_of_entry.resize(num_entries);
  std::fill(pos_of_entry.begin(), pos_of_entry.end(), -1);
  count_alleles.resize(num_entries);
  std::vector<std::vector<int>>::iterator it_count_alleles;
  for(it_count_alleles = count_alleles.begin(); it_count_alleles != count_alleles.end(); it_count_alleles++){
    (*it_count_alleles).resize(4);
  }

  //read first entry
  int ret = 1;
  if(chr == NULL){
		ret = sam_read1(fp_in, bamHdr, aln);
		chr = bamHdr->target_name[aln->core.tid];
	}
	if(strcmp(chr, contig.c_str()) != 0){
		while(strcmp(chr, contig.c_str()) != 0 && ret > 0){
			ret = sam_read1(fp_in, bamHdr, aln);
			chr = bamHdr->target_name[aln->core.tid];
		}
	}

  if(ret > 0){
    pos = aln->core.pos; //left most position of alignment in zero based coordianate (+1)
    chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
    if(contig != ""){
      if(strcmp(chr, contig.c_str()) != 0){
				std::cerr << chr << " " << contig << std::endl;
        std::cerr << "Error: contig names do not match" << std::endl;
        exit(1);
      }
    }else{
      contig = chr;
    }
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

}
