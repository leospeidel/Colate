#include <stdio.h>
#include <iostream>

#include "htslib.hpp"

#include "gzstream.hpp"
#include "data.hpp"
#include "sample.hpp"
#include "anc.hpp"
#include "mutations.hpp"
#include "cxxopts.hpp"

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include "coal_tree.hpp"
#include "coal_EM.hpp"

void
coal(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: num_bins, coal" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Calculating coalescence rates for (ancient) sample.." << std::endl;

  /////////////////////////////////
  //get TMRCA at each SNP

  ////////////////////////////////////////

  //decide on epochs
  int num_epochs = 0; 
  std::vector<double> epochs;

  std::ifstream is;
  std::string line;
  if(options.count("coal") > 0){
    is.open(options["coal"].as<std::string>());
    getline(is, line);
    getline(is, line);
    std::string tmp;
    for(int i = 0; i < line.size(); i++){
      if(line[i] == ' ' || line[i] == '\t'){
        epochs.push_back(stof(tmp));
        tmp = "";
        num_epochs++;
      }else{
        tmp += line[i];
      }
    }
    if(tmp != "") epochs.push_back(stof(tmp));

    assert(epochs[0] == 0);
    for(int e = 1; e < num_epochs; e++){
      std::cerr << epochs[e] << " ";
      assert(epochs[e] > epochs[e-1]);
    }

  }else{
    num_epochs = 30;
    if(options.count("num_bins") > 0){
      num_epochs = options["num_bins"].as<int>();
    }
    float years_per_gen = 28.0;
    if(options.count("years_per_gen")){
      years_per_gen = options["years_per_gen"].as<float>();
    }
    num_epochs++; 
    epochs.resize(num_epochs);

    epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }

  int num_bootstrap = 100;
  int block_size = 1000;

  coal_tree ct(epochs, num_bootstrap, block_size);

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file

  std::vector<std::string> chromosomes;
  if(options.count("chr") > 0){

    igzstream is_chr(options["chr"].as<std::string>());
    if(is_chr.fail()){
      std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
    }
    while(getline(is_chr, line)){
      chromosomes.push_back(line);
    }
    is_chr.close();

  }else{
    chromosomes.resize(1);
    chromosomes[0] = "1";
  }

  for(int chr = 0; chr < chromosomes.size(); chr++){

    std::cerr << "CHR " << chromosomes[chr] << ": ";
    AncMutIterators ancmut(options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".anc", options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".mut");
    float num_bases_tree_persists = 0.0;

    ct.update_ancmut(ancmut);

    int tree_count = 0, perc = -1;
    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    while(num_bases_tree_persists >= 0.0){
      if( (int) (((double)tree_count)/ancmut.NumTrees() * 100.0) > perc ){
        perc = (int) (((double)tree_count)/ancmut.NumTrees() * 100.0);
        std::cerr << "[" << perc << "%]\r";
      }
      tree_count++;
      ct.populate(mtr.tree);	
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    }

    std::cerr << std::endl;

  }

  ct.Dump(options["output"].as<std::string>() + ".coal");

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

int
test_vcf(cxxopts::Options&options){

  bool help = false;
  if(!options.count("input")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Parsing vcf file.." << std::endl;

  vcf_parser v(options["input"].as<std::string>());
  std::cerr << v.has_field("PL") << " " << v.has_field("foo") << std::endl;
  while(v.read_snp() == 0){
    v.extract_field("GT");
    for(int i = 0; i < v.n; i++){
      std::cerr << "(";
      for(int j = 0; j < v.num_entries_per_sample; j++){
        std::cerr << bcf_gt_allele(v.f[v.num_entries_per_sample*i+j]) << "|";
      }
      std::cerr << ") ";
    }
    std::cerr << std::endl;
  }

  return EXIT_SUCCESS;

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

int
test_bam(cxxopts::Options&options){

  bool help = false;
  if(!options.count("input")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Parsing bam file.." << std::endl;

  fasta ref_genome;
  ref_genome.Read(options["ref_genome"].as<std::string>());
  fasta anc_genome;
  //anc_genome.Read(options["anc_genome"].as<std::string>());

  bam_parser bam(options["target_bam"].as<std::string>(), options["ref_genome"].as<std::string>());

  haps mhaps((options["input"].as<std::string>() + ".haps.gz").c_str(), (options["input"].as<std::string>() + ".sample.gz").c_str());
  std::vector<char> seq(mhaps.GetN());
  int bp;

  for(int i = 0; i < 10; i++){

    mhaps.ReadSNP(seq, bp);

    if(bam.read_to_pos(bp)) break;

    std::cerr << bp << ": " << ref_genome.seq[bp-1] << " ";
    for(int i = 0; i < 4; i++){
      std::cerr << bam.count_alleles[(bp - 1) % bam.num_entries][i] << " ";
    }
    std::cerr << "| " << mhaps.ancestral << " " << mhaps.alternative << " " << seq[0] << " " << seq[1] << std::endl;

    if(0){
      int p;
      if(bam.pos_of_entry[p % bam.num_entries] == p){

        int num_alleles = 0, num_reads = 0;
        for(int i = 0; i < 4; i++){
          num_alleles += (bam.count_alleles[p % bam.num_entries][i] > 0);
          num_reads += bam.count_alleles[p % bam.num_entries][i];
        }
        if(num_alleles > 0){
          std::cerr << p << ": " << ref_genome.seq[p] << " ";
          for(int i = 0; i < 4; i++){
            std::cerr << bam.count_alleles[p % bam.num_entries][i] << " ";
          }
          std::cerr << " | " << num_alleles << " " << num_reads << std::endl;
        }
      }
    }

  }

  std::cerr << bam.pos << std::endl;
  std::cerr << bam.coverage/((double) bam.pos) << " " << bam.coverage_after_filter/((double) bam.pos) << std::endl;

  if(0){
    int entry_pos = 0;
    while(bam.read_entry() > 0){

      if(bam.pos - entry_pos > 5000){
        for(int p = entry_pos; p < bam.pos; p++){

          if(bam.pos_of_entry[p % bam.num_entries] == p){

            int num_alleles = 0;
            for(int i = 0; i < 4; i++){
              num_alleles += (bam.count_alleles[p % bam.num_entries][i] > 0);
            }
            if(num_alleles > 1){
              std::cerr << p << ": " << ref_genome.seq[p] << " ";
              for(int i = 0; i < 4; i++){
                std::cerr << bam.count_alleles[p % bam.num_entries][i] << " ";
              }
              std::cerr << " | " << num_alleles << std::endl;
            }
          }

        }
        entry_pos = bam.pos;
      }

    }
  }

  return EXIT_SUCCESS;

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

void
parse_vcf(std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_target_mask, double age, double C, std::mt19937& rng, std::vector<double>& age_shared_count, std::vector<double>& age_notshared_count){

  std::uniform_real_distribution<double> dist_unif(0,1);
  int num_age_bins = ((int) (log(1e8) * C))+1;
  int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0;
  bool target_flip = false;
  std::string ancestral, derived, ancestral_target, derived_target;
  Muts::iterator it_mut;

  bool has_mask = false;
  if(filename_target_mask.size() > 0) has_mask = true;

  for(int chr = 0; chr < filename_mut.size(); chr++){

    std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

    Mutations mut;
    mut.Read(filename_mut[chr]);
    vcf_parser target(filename_target[chr]);
    fasta mask;
    if(has_mask) mask.Read(filename_target_mask[chr]);

    bp_target = 0, bp_mut = 0;
    for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

      if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

        if((*it_mut).age_end >= age || (*it_mut).age_begin >= age){

          int i = 0;
          ancestral.clear();
          derived.clear();
          while((*it_mut).mutation_type[i] != '/' && i < (*it_mut).mutation_type.size()){
            ancestral.push_back((*it_mut).mutation_type[i]);
            i++;
          }
          i++;
          while(i < (*it_mut).mutation_type.size()){
            derived.push_back((*it_mut).mutation_type[i]);
            i++;
          }
          bp_mut = (*it_mut).pos;

          if(ancestral.size() > 0 && derived.size() > 0){

            //check if mutation exists in ref
            if(bp_target < bp_mut){
              while(target.read_snp() == 0){
                bp_target = target.rec->pos + 1;
                if(bp_target >= bp_mut) break;
              }
            }
            if(bp_target == bp_mut){

              ancestral_target = target.rec->d.allele[0];
              derived_target   = target.rec->d.allele[1];
              if( (ancestral_target == ancestral && derived_target == derived) || (ancestral_target == derived && derived_target == ancestral) ){

                target_flip = false;
                if( (ancestral_target == derived && derived_target == ancestral) ) target_flip = true;

                int bin_index;
                float num_samples = 10;
                double age_begin = (*it_mut).age_begin;

                //populate age_shared_count and age_notshared_count
                bool use = true;
                if(has_mask && bp_mut < mask.seq.size()){
                  if(mask.seq[bp_mut-1] == 'N') use = false;
                }

                if(use){

                  //possibility of sharing
                  target.extract_GT();

                  if(target.n*target.ploidy < 2){
                    std::cerr << "Need at least one diploid genome." << std::endl;
                    exit(1);
                  }

                  //N_ref       = target.n * target.ploidy/2;
                  //N_target    = N_ref;
                  N_ref = 0;
                  N_target = 0;

                  int DAF_target = 0;
                  int DAF_ref = 0;

                  bool biallelic = true;
                  int choose = (dist_unif(rng) < 0.5);
                  int k = 0;
                  for(int i = 0; i < target.n; i++){
                    for(int j = 0; j < target.ploidy; j++){
                      if(k % 2 == choose){
                        DAF_target += bcf_gt_allele(target.gt[target.ploidy*i+j]);
                        N_target++;
                      }else{
                        DAF_ref += bcf_gt_allele(target.gt[target.ploidy*i+j]);
                        N_ref++;
                      }
                      k++;
                      if(bcf_gt_allele(target.gt[target.ploidy*i+j]) > 1) biallelic = false;
                    }
                  }

                  if(biallelic){
                    assert(DAF_ref <= N_ref);
                    if(target_flip){
                      DAF_ref = N_ref - DAF_ref;
                      DAF_target = N_target - DAF_target;
                    }

                    if(DAF_ref > 0){

                      num_used_snps++;

                      bool skip = false;
                      //(*it_mut).age_end += age;
                      if((*it_mut).age_begin <= age){

                        if(1){
                          double age_begin2 = age_begin;
                          //age_begin2 = std::max(1.0*(*it_mut).age_begin, age);
                          age_begin2 = 0.0;
                          //age_begin2 = dist_unif(rng) * ((*it_mut).age_end - age_begin2) + age_begin2;
                          int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
                          int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
                          bin_index  = bin_index1 * num_age_bins + bin_index2;
                          age_shared_count[bin_index]    += DAF_target * DAF_ref;
                          age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref;

                          if(1){
                            int j = 0;
                            while(j < num_samples){
                              skip = false;
                              double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                              if(sampled_age < age) skip = true;
                              int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                              bin_index = bin_index_age * num_age_bins + bin_index_age;

                              if(!skip){
                                age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref/num_samples;
                              }
                              j++;
                            }
                          }
                        }

                      }else{

                        if(1){
                          int j = 0;
                          while(j < num_samples){
                            skip = false;
                            double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                            if(sampled_age < age) skip = true;
                            int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                            bin_index = bin_index_age * num_age_bins + bin_index_age;

                            if(!skip){
                              age_shared_count[bin_index]    += DAF_target * DAF_ref/num_samples;
                              age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref/num_samples;
                            }
                            j++;
                          }
                        }

                      }
                    }
                  }

                }
              }

            }
          }
        }
      }
    }

  }

  std::cerr << num_used_snps << std::endl << std::endl;

}

void
parse_vcfvcf(std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_ref, std::vector<std::string>& filename_target_mask, double age, double C, std::mt19937& rng, std::vector<double>& age_shared_count, std::vector<double>& age_notshared_count){

  std::uniform_real_distribution<double> dist_unif(0,1);
  int num_age_bins = ((int) (log(1e8) * C))+1;
  int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0;
  bool ref_flip = false, target_flip = false;
  std::string ancestral, derived, ancestral_ref, derived_ref, ancestral_target, derived_target;
  Muts::iterator it_mut;

  bool has_mask = false;
  if(filename_target_mask.size() > 0) has_mask = true;

  for(int chr = 0; chr < filename_mut.size(); chr++){

    std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

    Mutations mut;
    mut.Read(filename_mut[chr]);
    vcf_parser target(filename_target[chr]);
    vcf_parser ref(filename_ref[chr]);
    fasta mask;
    if(has_mask) mask.Read(filename_target_mask[chr]);

    bp_ref = 0, bp_target = 0, bp_mut = 0;
    for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

      if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

        if((*it_mut).age_end >= age || (*it_mut).age_begin >= age){

          int i = 0;
          ancestral.clear();
          derived.clear();
          while((*it_mut).mutation_type[i] != '/' && i < (*it_mut).mutation_type.size()){
            ancestral.push_back((*it_mut).mutation_type[i]);
            i++;
          }
          i++;
          while(i < (*it_mut).mutation_type.size()){
            derived.push_back((*it_mut).mutation_type[i]);
            i++;
          }
          bp_mut = (*it_mut).pos;

          if(ancestral.size() > 0 && derived.size() > 0){

            //check if mutation exists in ref
            if(bp_ref < bp_mut){
              while(ref.read_snp() == 0){
                bp_ref = ref.rec->pos + 1;
                if(bp_ref >= bp_mut) break;
              }
            }
            if(bp_ref == bp_mut){

              ancestral_ref = ref.rec->d.allele[0];
              derived_ref   = ref.rec->d.allele[1];
              if( (ancestral_ref == ancestral && derived_ref == derived) || (ancestral_ref == derived && derived_ref == ancestral) ){

                //assert(ancestral_ref == "A" || ancestral_ref == "C" || ancestral_ref == "G" || ancestral_ref == "T");
                //assert(derived_ref == "A" || derived_ref == "C" || derived_ref == "G" || derived_ref == "T");
                ref_flip = false;
                if( (ancestral_ref == derived && derived_ref == ancestral) ) ref_flip = true;

                //check if mutation exists in target
                if(bp_target < bp_mut){
                  while(target.read_snp() == 0){
                    bp_target = target.rec->pos + 1;
                    if(bp_target >= bp_mut) break;
                  }
                  target.extract_GT();
                  char** tmp = (target.rec -> d.allele);
                  int num_alleles = 0;
                  for(int i = 0; tmp[i] != NULL; i++){
                    num_alleles++;
                  }
                  N_target         = target.n * target.ploidy;
                  if(num_alleles > 0){
                    ancestral_target = target.rec->d.allele[0];
                  }else{
                    ancestral_target = "";
                  }
                  if(num_alleles > 1){
                    derived_target   = target.rec->d.allele[1];
                  }else{
                    derived_target   = "";
                  }
                }

                int bin_index;
                float num_samples = 10;
                //double age_begin = std::max(1.0*(*it_mut).age_begin, age);
                double age_begin = (*it_mut).age_begin;

                //populate age_shared_count and age_notshared_count
                bool use = true;
                if(has_mask && bp_mut < mask.seq.size()){
                  if(mask.seq[bp_mut-1] == 'N') use = false;
                }

                if(use){
                  if(bp_mut == bp_target){

                    //target could be fixed for reference
                    bool accept = false, fixed_for_ref = false;
                    target_flip = false;
                    if(ancestral_target != "" && derived_target == ""){
                      if(ancestral_target == ancestral || ancestral_target == derived) accept = true;
                      if(ancestral_target == derived) target_flip = true;
                      fixed_for_ref = true;
                    }else if( (ancestral_target == derived && derived_target == ancestral) || (ancestral_target == ancestral && derived_target == derived) ){
                      accept = true;
                      if( (ancestral_target == derived && derived_target == ancestral) ) target_flip = true;
                    }

                    int DAF_target = 0;
                    if( accept ){                  
                      bool biallelic = true;
                      for(int i = 0; i < target.n; i++){
                        for(int j = 0; j < target.ploidy; j++){
                          DAF_target += bcf_gt_allele(target.gt[target.ploidy*i+j]);
                          if(bcf_gt_allele(target.gt[target.ploidy*i+j]) > 1) biallelic = false;
                        }
                      }
                      if(!biallelic) accept = false;
                      if(fixed_for_ref && DAF_target != 0){
                        accept = false; //alt allele was not reported so this makes only sense if snp is fixed for ref in target
                      }
                    }

                    //possibility of sharing
                    if( accept ){

                      assert(DAF_target <= N_target);
                      ref.extract_GT();
                      N_ref       = ref.n * ref.ploidy;
                      int DAF_ref = 0;
                      bool biallelic = true;
                      for(int i = 0; i < ref.n; i++){
                        for(int j = 0; j < ref.ploidy; j++){
                          DAF_ref += bcf_gt_allele(ref.gt[ref.ploidy*i+j]);
                          if(bcf_gt_allele(ref.gt[ref.ploidy*i+j]) > 1) biallelic = false;
                        }
                      }

                      if(biallelic){
                        assert(DAF_ref <= N_ref);
                        assert(target_flip == ref_flip);
                        if(ref_flip){
                          DAF_ref = N_ref - DAF_ref;
                        }

                        if(DAF_ref > 0){

                          num_used_snps++;
                          if(target_flip){
                            DAF_target = N_target - DAF_target;
                          }
                          //assert(DAF_ref == (*it_mut).freq[0]);

                          bool skip = false;
                          //(*it_mut).age_end += age;
                          if((*it_mut).age_begin <= age){

                            if(1){
                              double age_begin2 = age_begin;
                              //age_begin2 = std::max(1.0*(*it_mut).age_begin, age);
                              age_begin2 = 0.0;
                              //age_begin2 = dist_unif(rng) * ((*it_mut).age_end - age_begin2) + age_begin2;
                              int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
                              int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
                              bin_index  = bin_index1 * num_age_bins + bin_index2;
                              age_shared_count[bin_index]    += DAF_target * DAF_ref;
                              //if(DAF_target > 0) age_shared_count[bin_index]++;
                              age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref;

                              if(1){
                                int j = 0;
                                //if(DAF_target > 0) age_begin = (*it_mut).age_end;
                                while(j < num_samples){
                                  skip = false;
                                  double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                                  if(sampled_age < age) skip = true;
                                  int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                                  bin_index = bin_index_age * num_age_bins + bin_index_age;

                                  if(!skip){
                                    //age_shared_count[bin_index]    += DAF_target * DAF_ref/num_samples;
                                    age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref/num_samples;
                                    //j++;
                                  }
                                  j++;
                                }
                              }
                            }

                          }else{

                            if(1){
                              int j = 0;
                              while(j < num_samples){
                                skip = false;
                                double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                                if(sampled_age < age) skip = true;
                                int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                                bin_index = bin_index_age * num_age_bins + bin_index_age;

                                if(!skip){
                                  age_shared_count[bin_index]    += DAF_target * DAF_ref/num_samples;
                                  age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref/num_samples;
                                  //j++;
                                }
                                j++;
                              }
                            }

                          }
                        }
                      }

                    }

                  }else{

                    ref.extract_GT();
                    N_ref       = ref.n * ref.ploidy;
                    int DAF_ref = 0, DAF_target = 0;
                    bool biallelic = true;
                    for(int i = 0; i < ref.n; i++){
                      for(int j = 0; j < ref.ploidy; j++){
                        DAF_ref += bcf_gt_allele(ref.gt[ref.ploidy*i+j]);
                        if(bcf_gt_allele(ref.gt[ref.ploidy*i+j]) > 1) biallelic = false;
                      }
                    }
                    assert(DAF_ref <= N_ref);

                    //assume that target is ref
                    DAF_target = 0;
                    if(ref_flip){
                      DAF_ref = N_ref - DAF_ref;
                      DAF_target = N_target - DAF_target;
                    }

                    if(biallelic && DAF_ref > 0){
                      num_used_snps++;

                      bool skip = false;
                      //(*it_mut).age_end += age;
                      if((*it_mut).age_begin <= age){

                        if(1){
                          double age_begin2 = age_begin;
                          age_begin2 = std::max(1.0*(*it_mut).age_begin, age);
                          age_begin2 = 0.0;
                          //age_begin2 = dist_unif(rng) * ((*it_mut).age_end - age_begin2) + age_begin2;
                          int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
                          int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
                          bin_index  = bin_index1 * num_age_bins + bin_index2;
                          age_shared_count[bin_index]    += DAF_target * DAF_ref;
                          //if(DAF_target > 0) age_shared_count[bin_index]++;
                          age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref;

                          if(1){
                            int j = 0;
                            //if(DAF_target > 0) age_begin = (*it_mut).age_end;
                            while(j < num_samples){
                              skip = false;
                              double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                              if(sampled_age < age) skip = true;
                              int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                              bin_index = bin_index_age * num_age_bins + bin_index_age;

                              if(!skip){
                                //age_shared_count[bin_index]    += DAF_target * DAF_ref/num_samples;
                                age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref/num_samples;
                                //j++;
                              }
                              j++;
                            }
                          }
                        }

                      }else{

                        if(1){
                          int j = 0;
                          while(j < num_samples){
                            skip = false;
                            double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                            if(sampled_age < 1.1*age) skip = true;
                            int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                            bin_index = bin_index_age * num_age_bins + bin_index_age;

                            if(!skip){
                              age_shared_count[bin_index]    += DAF_target * DAF_ref/num_samples;
                              age_notshared_count[bin_index] += (N_target - DAF_target) * DAF_ref/num_samples;
                              //j++;
                            }
                            j++;
                          }
                        }

                      }
                    }

                  }

                }
              }

            }
          }
        }
      }
    }

  }

  std::cerr << num_used_snps << std::endl << std::endl;

}

void
parse_bamvcf(std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_ref, std::vector<std::string>& filename_ref_genome, double age, double C, std::mt19937& rng, std::vector<double>& age_shared_count, std::vector<double>& age_notshared_count){

  std::uniform_real_distribution<double> dist_unif(0,1);
  int num_age_bins = ((int) (log(1e8) * C))+1;
  int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0, num_used_snps2 = 0;
  bool ref_flip = false, target_flip = false;
  std::string ancestral, derived, ancestral_ref, derived_ref, ancestral_target, derived_target;
  Muts::iterator it_mut;
  for(int chr = 0; chr < filename_mut.size(); chr++){

    std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

    Mutations mut;
    mut.Read(filename_mut[chr]);
    bam_parser target(filename_target[chr], filename_ref_genome[chr]);
    vcf_parser ref(filename_ref[chr]);	

    bp_ref = 0, bp_target = 0, bp_mut = 0;
    for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

      if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

        if((*it_mut).age_end >= 2*age){

          int i = 0;
          ancestral.clear();
          derived.clear();
          while((*it_mut).mutation_type[i] != '/' && i < (*it_mut).mutation_type.size()){
            ancestral.push_back((*it_mut).mutation_type[i]);
            i++;
          }
          i++;
          while(i < (*it_mut).mutation_type.size()){
            derived.push_back((*it_mut).mutation_type[i]);
            i++;
          }
          bp_mut = (*it_mut).pos;

          if(ancestral.size() > 0 && derived.size() > 0){

            //check if mutation exists in ref
            if(bp_ref < bp_mut){
              while(ref.read_snp() == 0){
                bp_ref = ref.rec->pos + 1;
                if(bp_ref >= bp_mut) break;
              }
            }
            if(bp_ref == bp_mut){

              num_used_snps2++;
              ancestral_ref = ref.rec->d.allele[0];
              derived_ref   = ref.rec->d.allele[1];
              if( (ancestral_ref == ancestral && derived_ref == derived) || (ancestral_ref == derived && derived_ref == ancestral) ){

                ref_flip = false;
                if( (ancestral_ref == derived && derived_ref == ancestral) ) ref_flip = true;

                //check if mutation exists in target     
                int bp_target = bp_mut - 1;
                target.read_to_pos(bp_target);
                int num_alleles = 0;
                bool use_snp = false;
                int DAF_target = 0, AAF_target = 0;
                if(target.pos_of_entry[bp_target % target.num_entries] == bp_target){
                  int num_reads = 0;
                  for(int i = 0; i < 4; i++){
                    num_reads += target.count_alleles[bp_target % target.num_entries][i];
                    num_alleles += (target.count_alleles[bp_target % target.num_entries][i] > 0);
                  }

                  if(ancestral == "A"){
                    AAF_target = target.count_alleles[bp_target % target.num_entries][0];
                  }else if(ancestral == "C"){
                    AAF_target = target.count_alleles[bp_target % target.num_entries][1];
                  }else if(ancestral == "G"){
                    AAF_target = target.count_alleles[bp_target % target.num_entries][2];
                  }else if(ancestral == "T"){
                    AAF_target = target.count_alleles[bp_target % target.num_entries][3];
                  }
                  if(derived == "A"){
                    DAF_target = target.count_alleles[bp_target % target.num_entries][0];
                  }else if(derived == "C"){
                    DAF_target = target.count_alleles[bp_target % target.num_entries][1];
                  }else if(derived == "G"){
                    DAF_target = target.count_alleles[bp_target % target.num_entries][2];
                  }else if(derived == "T"){
                    DAF_target = target.count_alleles[bp_target % target.num_entries][3];
                  }

                  if( (AAF_target > 0 || DAF_target > 0)){
                    if(num_alleles <= 2) use_snp = true;
                    //DAF_target /= (DAF_target + AAF_target);
                    //AAF_target /= (DAF_target + AAF_target);
                    //std::cerr << bp_mut << " " << ancestral << " " << derived << " " << (*it_mut).age_begin << " " << (*it_mut).age_end << " | " << num_reads << " " << AAF_target << " " << DAF_target  << std::endl;
                  }

                }

                int bin_index;
                float num_samples = 10;
                //double age_begin = std::max(1.0*(*it_mut).age_begin, age);
                double age_begin = (*it_mut).age_begin;

                //populate age_shared_count and age_notshared_count
                if(use_snp){

                  //possibility of sharing
                  ref.extract_GT();
                  N_ref       = ref.n * ref.ploidy;
                  int DAF_ref = 0;
                  for(int i = 0; i < ref.n; i++){
                    for(int j = 0; j < ref.ploidy; j++){
                      DAF_ref += bcf_gt_allele(ref.gt[ref.ploidy*i+j]);
                    }
                  }
                  assert(DAF_ref <= N_ref);
                  if(ref_flip){
                    DAF_ref = N_ref - DAF_ref;
                  }

                  if(DAF_ref > 0){

                    //DAF_ref /= N_ref;
                    /*
                       int a = DAF_target, b = AAF_target;
                       DAF_target /= (DAF_target + AAF_target);
                       AAF_target /= (DAF_target + AAF_target);
                       */
                    /*
                       if(dist_unif(rng) < DAF_target){
                       DAF_target = 1;
                       AAF_target = 0;
                       }else{
                       DAF_target = 0;
                       AAF_target = 1;
                       }
                       */
                    /*
                       if(DAF_target > 0.66){
                       DAF_target = 2;
                       AAF_target = 0;
                       }else if(DAF_target > 0.33){
                       DAF_target = 1;
                       AAF_target = 1;
                       }else{
                       DAF_target = 0;
                       AAF_target = 2;
                       }
                       std::cerr << bp_mut << " " << ancestral << " " << derived << " " << a << " " << b << " " << DAF_target << std::endl;
                       */


                    //std::cerr << AAF_target << " " << DAF_target << " " << target.count_alleles[bp_mut % target.num_entries][0]  << " " << target.count_alleles[bp_mut % target.num_entries][1] << " " << target.count_alleles[bp_mut % target.num_entries][2] << " " << target.count_alleles[bp_mut % target.num_entries][3] << std::endl;

                    bool skip = false;
                    //(*it_mut).age_end += age;
                    if((*it_mut).age_begin <= age){

                      if(1){
                        double age_begin2 = age_begin;
                        age_begin2 = 0.0;
                        age_begin2 = std::max(1.0*(*it_mut).age_begin, age);
                        int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
                        int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
                        bin_index  = bin_index1 * num_age_bins + bin_index2;
                        age_shared_count[bin_index]    += DAF_target * DAF_ref;
                        //age_notshared_count[bin_index] += AAF_target * DAF_ref;

                        if(1){
                          int j = 0;
                          while(j < num_samples){
                            skip = false;
                            double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                            if(sampled_age < age && DAF_target > 0) skip = true;
                            int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                            bin_index = bin_index_age * num_age_bins + bin_index_age;

                            if(!skip){
                              age_notshared_count[bin_index] += AAF_target * DAF_ref/num_samples;
                              j++;
                            }
                          }
                        }
                      }

                    }else{

                      int j = 0;
                      while(j < num_samples){
                        skip = false;
                        double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
                        if(sampled_age < age && DAF_target > 0) skip = true;
                        int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
                        bin_index = bin_index_age * num_age_bins + bin_index_age;

                        if(!skip){
                          age_shared_count[bin_index]    += DAF_target * DAF_ref/num_samples;
                          age_notshared_count[bin_index] += AAF_target * DAF_ref/num_samples;
                          j++;
                        }
                      }

                    }

                    //num_used_snps++;
                  }

                }

              }

            }

          }
        }
      }
    }
    std::cerr << "Coverage: " << ((double) target.coverage)/target.ref_genome.seq.size() << " " << ((double) target.coverage_after_filter)/target.ref_genome.seq.size() << std::endl; 
  }

}


void
mut(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("mut") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mut, target_vcf, reference_vcf, bins, output. Optional: target_bam, ref_genome, target_age, reference_age, coal" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Calculating coalescence rates for (ancient) sample.." << std::endl;

  int regularise = 2;
  if(options.count("regularise") > 0){ 
    regularise = options["regularise"].as<int>();
  }

  bool is_ancient = false;
  double age;
  double log_10 = std::log(10);

  double target_age = 0;
  double ref_age    = 0;
  if(options.count("target_age") > 0) target_age = std::stof(options["target_age"].as<std::string>());
  if(options.count("ref_age") > 0)    ref_age    = std::stof(options["reference_age"].as<std::string>());
  assert(target_age >= 0.0);
  assert(ref_age >= 0.0);

  double years_per_gen = 28.0;
  if(options.count("years_per_gen") > 0){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  age = std::max(target_age, ref_age)/years_per_gen;
  std::cerr << age << std::endl; 
  if(age > 0.0) is_ancient = true;

  std::string line;
  ////////////////////////////////////////

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 10;
  int num_age_bins = ((int) (log(1e8) * C))+1;
  std::vector<double> age_bin(num_age_bins, 0.0);
  std::vector<double>::iterator it_age_bin = age_bin.begin();
  int bin = 0;
  *it_age_bin = 0.0;
  it_age_bin++;
  for(; bin < num_age_bins-1; bin++){
    *it_age_bin =  exp(bin/C)/10.0;
    it_age_bin++;
  }
  std::vector<double> age_shared_count(num_age_bins*num_age_bins, 0), age_notshared_count(num_age_bins*num_age_bins, 0);

  ///////////

  std::mt19937 rng;
  int seed = std::time(0) + getpid();
  if(options.count("seed") > 0){
    seed = options["seed"].as<int>();
  }
  rng.seed(seed);

  igzstream is_mat(options["output"].as<std::string>() + ".colate_mat");
  if(is_mat.fail()){

    ///////////
    std::vector<std::string> filename_mut, filename_target, filename_ref, filename_ref_genome, filename_target_mask;

    if(options.count("target_vcf") && options.count("reference_vcf")){
      if(options.count("chr") > 0){

        igzstream is_chr(options["chr"].as<std::string>());
        if(is_chr.fail()){
          std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
        }
        while(getline(is_chr, line)){
          filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
          filename_target.push_back(options["target_vcf"].as<std::string>() + "_chr" + line + ".bcf");
          filename_ref.push_back(options["reference_vcf"].as<std::string>() + "_chr" + line + ".bcf");
          if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
        }
        is_chr.close();

      }else{
        filename_mut.push_back(options["mut"].as<std::string>());
        filename_target.push_back(options["target_vcf"].as<std::string>());
        filename_ref.push_back(options["reference_vcf"].as<std::string>());
        if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>());
      }
      parse_vcfvcf(filename_mut, filename_target, filename_ref, filename_target_mask, age, C, rng, age_shared_count, age_notshared_count);

    }else if(options.count("target_vcf")){
      if(options.count("chr") > 0){

        igzstream is_chr(options["chr"].as<std::string>());
        if(is_chr.fail()){
          std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
        }
        while(getline(is_chr, line)){
          filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
          filename_target.push_back(options["target_vcf"].as<std::string>() + "_chr" + line + ".bcf");
          if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
        }
        is_chr.close();

      }else{
        filename_mut.push_back(options["mut"].as<std::string>());
        filename_target.push_back(options["target_vcf"].as<std::string>());
        if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>());
      }
      parse_vcf(filename_mut, filename_target, filename_target_mask, age, C, rng, age_shared_count, age_notshared_count);

    }else if(options.count("target_bam") && options.count("reference_vcf")){
      if(options.count("chr") > 0){

        igzstream is_chr(options["chr"].as<std::string>());
        if(is_chr.fail()){
          std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
        }
        while(getline(is_chr, line)){
          filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
          filename_target.push_back(options["target_bam"].as<std::string>() + "_chr" + line + ".bam");
          filename_ref.push_back(options["reference_vcf"].as<std::string>() + "_chr" + line + ".bcf");
          filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
        }
        is_chr.close();

      }else{
        filename_mut.push_back(options["mut"].as<std::string>());
        filename_target.push_back(options["target_bam"].as<std::string>());
        filename_ref.push_back(options["reference_vcf"].as<std::string>());
        filename_ref_genome.push_back(options["ref_genome"].as<std::string>());
      }
      parse_bamvcf(filename_mut, filename_target, filename_ref, filename_ref_genome, age, C, rng, age_shared_count, age_notshared_count);

    }

    std::vector<double> F(num_age_bins, 0.0), G(num_age_bins, 0.0);
    double fcount = 0.0, gcount = 0.0;
    int bin = 0;
    while(age_bin[bin] <= age) bin++;
    //std::cerr << bin << " " << age_bin[bin]*28 << std::endl;
    for(; bin < num_age_bins; bin++){
      //std::cerr << age_shared_count[bin] << std::endl;
      if(age_shared_count[bin] > 0){
        F[bin] = F[bin-1] + age_shared_count[bin]/(age_shared_count[bin] + age_notshared_count[bin]) * age_bin[bin];
        //F[bin] = F[bin-1] + age_shared_count[bin];
      }else{
        F[bin] = F[bin-1];
      }
      if(age_notshared_count[bin] > 0){
        G[bin] = G[bin-1] + age_notshared_count[bin]/(age_shared_count[bin] + age_notshared_count[bin]) * age_bin[bin];
        //G[bin] = G[bin-1] + age_notshared_count[bin];
      }else{
        G[bin] = G[bin-1];
      }
      fcount += age_shared_count[bin];
      gcount += age_notshared_count[bin];
      age_shared_count[bin] = 0.0;
      age_notshared_count[bin] = 0.0;
    }
    if(1){
      if(F[num_age_bins-1] > 0){
        for(bin = 1; bin < num_age_bins; bin++){
          F[bin]  = F[bin]/F[num_age_bins-1];
          F[bin] *= fcount;
          int bin_index = bin*num_age_bins+bin;
          age_shared_count[bin_index] += F[bin]-F[bin-1];
        }
      }
    }
    if(0){
      if(G[num_age_bins-1] > 0){
        for(bin = 1; bin < num_age_bins-1; bin++){
          G[bin] = G[bin]/G[num_age_bins-1];
          G[bin] *= gcount;
          int bin_index = bin*num_age_bins+bin;
          age_notshared_count[bin_index] += G[bin]-G[bin-1];
        }
      }
    }

    //print matrix
    std::ofstream os_mat(options["output"].as<std::string>() + ".colate_mat");  
    for(int bin = 0; bin < num_age_bins; bin++){
      os_mat << age_bin[bin] << " "; 
    }
    os_mat << "\n";
    for(int bin1 = 0; bin1 < num_age_bins; bin1++){
      for(int bin2 = 0; bin2 < num_age_bins; bin2++){
        int bin = bin1*num_age_bins + bin2;
        os_mat << age_shared_count[bin] << " ";
      }
      os_mat << "\n";
    }
    for(int bin1 = 0; bin1 < num_age_bins; bin1++){
      for(int bin2 = 0; bin2 < num_age_bins; bin2++){
        int bin = bin1*num_age_bins + bin2;
        os_mat << age_notshared_count[bin] << " ";
      }
      os_mat << "\n";
    }
    os_mat.close();

  }else{
    std::cerr << "Loading precomputed file " << options["output"].as<std::string>() << ".colate_mat" << std::endl;
    for(int bin = 0; bin < num_age_bins; bin++){
      is_mat >> age_bin[bin]; 
    }
    for(int bin1 = 0; bin1 < num_age_bins; bin1++){
      for(int bin2 = 0; bin2 < num_age_bins; bin2++){
        int bin = bin1*num_age_bins + bin2;
        is_mat >> age_shared_count[bin];
      }
    }
    for(int bin1 = 0; bin1 < num_age_bins; bin1++){
      for(int bin2 = 0; bin2 < num_age_bins; bin2++){
        int bin = bin1*num_age_bins + bin2;
        is_mat >> age_notshared_count[bin];
      }
    }
  }

  ///////////////////////////////////////

  //decide on epochs
  int num_epochs = 0; 
  std::vector<double> epochs;
  std::ifstream is;
  if(options.count("coal") > 0){
    is.open(options["coal"].as<std::string>());
    getline(is, line);
    getline(is, line);
    std::string tmp;
    for(int i = 0; i < line.size(); i++){
      if(line[i] == ' ' || line[i] == '\t'){
        epochs.push_back(stof(tmp));
        tmp = "";
        num_epochs++;
      }else{
        tmp += line[i];
      }
    }
    if(tmp != "") epochs.push_back(stof(tmp));

    assert(epochs[0] == 0);
    for(int e = 1; e < num_epochs; e++){
      std::cerr << epochs[e] << " ";
      assert(epochs[e] > epochs[e-1]);
    }

  }else{ 

    double min_epoch = 0.0;
    double count = 0.0;
    int bin = 0;
    for(; bin < num_age_bins; bin++){
      count += age_shared_count[bin] + age_shared_count[bin*num_age_bins+bin];
      if(count > 100) break;
    }
    min_epoch = (((int)(log10(age_bin[bin]*years_per_gen)*100)))/100.0;
    double log_age = std::log(age * years_per_gen)/log_10;

    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    //epoch_lower = std::max(min_epoch, (double) std::stof(tmp));
    //std::cerr << epoch_lower << std::endl;
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      ep++;
    }
    epoch_boundary = epoch_lower;
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0){
          epochs.push_back(age);
        }
        if(std::fabs(log_age - epoch_boundary) > 1e-3){
          epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        }
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10*epochs[epochs.size()-1])/years_per_gen );
    num_epochs = epochs.size();	

  }

  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)
  double initial_coal_rate = 1.0/20000.0;
  std::vector<double> coal_rates(num_epochs, initial_coal_rate), coal_rates_num(num_epochs, 0.0), coal_rates_denom(num_epochs, 0.0);
  if(options.count("coal") > 0){
    double dummy;
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      is >> coal_rates[e];
      std::cerr << coal_rates[e] << " "; 
    }
    std::cerr << std::endl;
  }


  ////////////////////////////////////////

  std::cerr << "Maximising likelihood using EM.. " << std::endl;

  std::ofstream os_log(options["output"].as<std::string>() + ".log");

  //coal_EM EM(epochs, coal_rates);
  double gamma = 1;
  std::vector<double> gradient(num_epochs, 0.0), gradient_prev(num_epochs, 0.0), coal_rates_prev(num_epochs, 0.0), coal_rates_nesterov(num_epochs, 0.0);
  coal_rates_nesterov = coal_rates;

  int max_iter = 10000;
  int perc = -1;
  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      std::cerr << "[" << perc << "%]\r";
    }

    //EM.UpdateCoal(coal_rates_nesterov);
    if(is_ancient){
      coal_rates_nesterov[0] = 0;
    }
    coal_EM EM(epochs, coal_rates_nesterov);
    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;

    double count = 0;
    std::vector<double> num(num_epochs,0.0), denom(num_epochs,0.0);
    for(int bin1 = 0; bin1 < num_age_bins; bin1++){
      for(int bin2 = bin1; bin2 < num_age_bins; bin2++){

        int bin = bin1*num_age_bins + bin2;
        if(age_shared_count[bin] > 0){ 
          count = age_shared_count[bin];
          double logl = EM.EM_shared(age_bin[bin1], age_bin[bin2], num, denom);
          log_likelihood += count * logl;
          for(int e = 0; e < num_epochs; e++){
            assert(!std::isnan(num[e]));
            assert(!std::isnan(denom[e]));
            assert(num[e] >= 0.0);
            assert(denom[e] >= 0.0);
            coal_rates_num[e]   += count * num[e];
            coal_rates_denom[e] += count * denom[e];
          }
        }
        if(age_notshared_count[bin] > 0){ 
          count = age_notshared_count[bin];
          double logl = EM.EM_notshared(age_bin[bin1], age_bin[bin2], num, denom);
          log_likelihood += count * logl;
          for(int e = 0; e < num_epochs; e++){
            assert(!std::isnan(num[e]));
            assert(!std::isnan(denom[e]));
            assert(num[e] >= 0.0);
            assert(denom[e] >= 0.0);
            coal_rates_num[e]   += count * num[e];
            coal_rates_denom[e] += count * denom[e];
          }
        }

      }	
    }

    bool is_EM = (regularise == 2);
    bool is_regular = (regularise == 1);

    if(is_regular){
      double alpha = 1.0;
      double beta  = 4e-4; 
      coal_rates_denom[0] += alpha*beta * tanh( (beta/coal_rates[1] - beta/coal_rates[0]) )/(coal_rates[0] * coal_rates[0]);
      for(int e = 1; e < num_epochs-1; e++){
        coal_rates_denom[e] += alpha*beta * ( tanh( (beta/coal_rates[e+1] - beta/coal_rates[e]) ) + tanh( (beta/coal_rates[e-1] - beta/coal_rates[e]) ) )/(coal_rates[e] * coal_rates[e]);
      }
      coal_rates_denom[num_epochs-1] += alpha*beta*tanh((beta/coal_rates[num_epochs-2] - beta/coal_rates[num_epochs-1]))/(coal_rates[num_epochs-1] * coal_rates[num_epochs-1]);
    }

    //calculate the gradient
    if(!is_EM){
      if(iter == 0){
        for(int e = 0; e < num_epochs; e++){
          if(coal_rates_nesterov[e] != 0.0 && coal_rates_num[e] != 0.0 && coal_rates_denom[e] != 0.0){
            gradient[e]        = (coal_rates_num[e]/coal_rates_nesterov[e] - coal_rates_denom[e]);
          }else{
            gradient[e]        = 0.0;
          }
        }
      }else{
        for(int e = 0; e < num_epochs; e++){
          if(coal_rates_nesterov[e] != 0.0 && coal_rates_num[e] != 0.0 && coal_rates_denom[e] != 0.0){
            gradient[e]        = (coal_rates_num[e]/coal_rates_nesterov[e] - coal_rates_denom[e]);
          }else{
            gradient[e]        = 0.0;
          }
        }
      }
      gamma = 1e-15;
    }

    //update coal_rates
    for(int e = 0; e < num_epochs; e++){

      if(!is_EM){
        coal_rates_prev[e] = coal_rates[e];
      }

      if(coal_rates_num[e] == 0){
        if(e > 0){
          coal_rates[e]          = coal_rates[e-1];
          coal_rates_nesterov[e] = coal_rates_nesterov[e-1];
        }else{
          coal_rates[e] = 0;
          coal_rates_nesterov[e] = 0;
        }
      }else if(coal_rates_denom[e] == 0){
      }else{

        if(is_EM){
          coal_rates[e] = coal_rates_num[e]/coal_rates_denom[e];
          coal_rates_nesterov[e] = coal_rates[e];
        }else{
          coal_rates[e]  += 0.9*gradient_prev[e] + gamma * gradient[e];
          coal_rates_nesterov[e] = coal_rates[e] + 0.9*(0.9*gradient_prev[e] + gamma * gradient[e]);
          //coal_rates[e] += gamma * gradient[e]; 
          //coal_rates_nesterov[e] = coal_rates[e];
        }
        if(1){
          if(coal_rates[e] < 1e-9){
            coal_rates[e] = 1e-9;
          }
          if(coal_rates_nesterov[e] < 1e-9){
            coal_rates_nesterov[e] = 1e-9;
          }
        }
        assert(coal_rates[e] >= 0.0);
        assert(coal_rates_nesterov[e] >= 0.0);

      }

      if(!is_EM){
        //gradient_prev[e]   = gradient[e];
        gradient_prev[e] = 0.9*gradient_prev[e] + gamma * gradient[e];
      }

    }

    os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << std::endl;
    std::fill(coal_rates_num.begin(), coal_rates_num.end(), 0.0);
    std::fill(coal_rates_denom.begin(), coal_rates_denom.end(), 0.0);

  }
  os_log.close();

  std::ofstream os(options["output"].as<std::string>() + ".coal");

  os << "0\n";
  for(int e = 0; e < num_epochs; e++){
    os << epochs[e] << " ";
  }
  os << std::endl;

  os << "0 " << "0" << " ";
  for(int e = 0; e < num_epochs; e++){
    os << coal_rates[e] << " ";
  }
  os << std::endl;

  os.close();

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

void
make_mut_incl_out(cxxopts::Options& options){

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, haps, sample, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Make mut file including fixed mutations." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Adding fixed mutations to mut file.." << std::endl;

  double outgroup_age = 10e6/28;

  haps mhaps(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  int N = mhaps.GetN();
  int L = mhaps.GetL();
  std::vector<char> sequence(N);
  int bp = -1, DAF;

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
  float num_bases_snp_persists = 0.0;
  std::vector<float> coords;
  double tmrca = 0.0;
  int root = 2*ancmut.NumTips()-2;

  //Mutations mut_out;
  //mut_out.Read(options["input"].as<std::string>() + ".mut");

  Mutations mut_combined;
  mut_combined.info.resize(L);

  int L_ref = ancmut.NumSnps();
  //int L_out = mut_out.info.size();

  std::string ancestral, derived;
  int snp_ref = 0, snp_out = 0, snp_hap = 0, snp_comb = 0;
  num_bases_snp_persists = ancmut.FirstSNP(mtr, it_mut);
  mtr.tree.GetCoordinates(coords);
  tmrca = coords[root];

  int tree_count = (*it_mut).tree;
  for(; snp_hap < L; snp_hap++){

    mhaps.ReadSNP(sequence, bp);
    DAF = 0;
    for(std::vector<char>::iterator it_seq = sequence.begin(); it_seq != sequence.end(); it_seq++){
      DAF += (*it_seq == '1');
    }

    //check if snp exists in mut, and whether it has freq > 0 in ref
    if(snp_ref < L_ref){
      while((*it_mut).pos < bp){
        num_bases_snp_persists = ancmut.NextSNP(mtr, it_mut);
        snp_ref++;
        if(snp_ref == L_ref) break;
      }
    }

    if(tree_count < (*it_mut).tree){
      tree_count = (*it_mut).tree;
      mtr.tree.GetCoordinates(coords);
      tmrca = coords[root];
    }
    //std::cerr << tmrca << " " << coords[root] << " " << root << std::endl;
    //if(snp_out < L_out){
    //	while(mut_out.info[snp_out].pos < bp){
    //		snp_out++;
    //		if(snp_out == L_out) break;
    //	}
    //}

    if((*it_mut).pos == bp && DAF > 0 && DAF < N){

      if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1){

        std::string ancestral, derived, haps_ancestral = mhaps.ancestral, haps_derived = mhaps.alternative;
        int i = 0;
        ancestral.clear();
        derived.clear();
        while((*it_mut).mutation_type[i] != '/' && i < (*it_mut).mutation_type.size()){
          ancestral.push_back((*it_mut).mutation_type[i]);
          i++;
        }
        i++;
        while(i < (*it_mut).mutation_type.size()){
          derived.push_back((*it_mut).mutation_type[i]);
          i++;
        }

        if( (ancestral == haps_ancestral && derived == haps_derived) || (derived == haps_ancestral && ancestral == haps_derived) ){

          if((derived == haps_ancestral && ancestral == haps_derived)) DAF = N-DAF;

          if((*it_mut).age_end > 0){

            //if segregating, copy over
            mut_combined.info[snp_comb] = (*it_mut);
            //if(mut_combined.info[snp_comb].branch.size() == 1){
            //  assert( *mut_combined.info[snp_comb].branch.begin() <= root );
            //}
            mut_combined.info[snp_comb].tree = tree_count;
            mut_combined.info[snp_comb].pos = bp;
            mut_combined.info[snp_comb].freq.clear();
            mut_combined.info[snp_comb].freq.push_back(DAF);
            if(snp_comb > 0) mut_combined.info[snp_comb-1].dist = bp - mut_combined.info[snp_comb-1].pos; 
            mut_combined.info[snp_comb].snp_id = snp_comb;

            //if(DAF == 1) assert(mut_combined.info[snp_comb].age_begin == 0.0);
            if(mut_combined.info[snp_comb].age_begin == 0.0){
              assert(DAF == 1);
            }

            snp_comb++;
          }

        }

      }

    }else{
      //otherwise copy from trees with outgroup
      if(0){
        /*
           if(DAF == N && bp == mut_out.info[snp_out].pos){
           if(tmrca <= mut_out.info[snp_out].age_end){
           mut_combined.info[snp_hap].branch.resize(1);
           mut_combined.info[snp_hap].branch[0] = root;
           mut_combined.info[snp_hap].age_begin = tmrca;
           mut_combined.info[snp_hap].age_end   = mut_out.info[snp_out].age_end;
           assert(mut_combined.info[snp_hap].age_begin <= mut_combined.info[snp_hap].age_end);
           }
           }else{
           mut_combined.info[snp_hap].age_begin = 0;
           mut_combined.info[snp_hap].age_end = 0;
           mut_combined.info[snp_hap].branch.clear();
           }
           */
      }else{     

        if(DAF == N){
          if(tmrca <= outgroup_age){
            mut_combined.info[snp_comb].branch.resize(1);
            mut_combined.info[snp_comb].branch[0] = root;
            mut_combined.info[snp_comb].age_begin = tmrca;
            mut_combined.info[snp_comb].age_end   = outgroup_age;
            assert(mut_combined.info[snp_comb].age_begin <= mut_combined.info[snp_comb].age_end);

            mut_combined.info[snp_comb].tree = tree_count;
            mut_combined.info[snp_comb].pos = bp;
            mut_combined.info[snp_comb].freq.push_back(DAF);

            ancestral = mhaps.ancestral;
            derived   = mhaps.alternative;
            mut_combined.info[snp_comb].mutation_type = ancestral + "/" + derived;
            if(snp_comb > 0) mut_combined.info[snp_comb-1].dist = bp - mut_combined.info[snp_comb-1].pos; 
            mut_combined.info[snp_comb].snp_id = snp_comb;

            snp_comb++;
          }
        }
      }
    }

  }

  mut_combined.info.resize(snp_comb);
  mut_combined.Dump(options["output"].as<std::string>());

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

void
preprocess_mut(cxxopts::Options& options){

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("reference_vcf") || !options.count("ref_genome") || !options.count("anc_genome") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, reference_vcf, ref_genome, anc_genome, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Make mut file including fixed mutations." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Adding fixed mutations to mut file.." << std::endl;

  double outgroup_age = 10e6/28;

  fasta ref_genome, anc_genome, mask;
  ref_genome.Read(options["ref_genome"].as<std::string>());
  anc_genome.Read(options["anc_genome"].as<std::string>());
  mask.Read(options["mask"].as<std::string>());

  vcf_parser vcf(options["reference_vcf"].as<std::string>());

  int bp = -1, bp_prev, DAF;

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
  float num_bases_snp_persists = 0.0;
  std::vector<float> coords;
  double tmrca = 0.0;
  int root = 2*ancmut.NumTips()-2;

  int L_ref = ancmut.NumSnps();
  Mutations mut_combined;
  mut_combined.info.resize(10*L_ref);

  std::string ancestral, derived;
  int snp_ref = 0, snp_out = 0, snp_hap = 0, snp_comb = 0;
  num_bases_snp_persists = ancmut.FirstSNP(mtr, it_mut);
  mtr.tree.GetCoordinates(coords);
  tmrca = coords[root];

  bool biallelic;
  int tree_count = (*it_mut).tree;
  while(vcf.read_snp() == 0){

    if(mut_combined.info.size() - snp_comb < L_ref) mut_combined.info.resize(mut_combined.info.size() + L_ref);

    vcf.extract_GT();
    int N = vcf.ploidy*vcf.n;
    if(N != ancmut.NumTips()){
      std::cerr << N << " " << ancmut.NumTips() << std::endl;
      std::cerr << "Number of samples in vcf and anc/mut don't match" << std::endl;
      exit(1);
    }

    biallelic = true;
    DAF = 0;
    for(int i = 0; i < vcf.n; i++){
      for(int j = 0; j < vcf.ploidy; j++){
        DAF += bcf_gt_allele(vcf.gt[vcf.ploidy*i+j]);
        if(bcf_gt_allele(vcf.gt[vcf.ploidy*i+j]) > 1) biallelic = false;
      }
    }
    bp_prev = bp;
    bp = vcf.rec -> pos;

    //go through all positions between bp_prev and bp and check if ref == anc 
    for(int bp_tmp = bp_prev+1; bp_tmp < bp; bp_tmp++){
      if(bp_tmp < mask.seq.size() && bp_tmp < anc_genome.seq.size() && bp_tmp < ref_genome.seq.size()){
        if(mask.seq[bp_tmp] == 'P'){
          //assume all are equal to ref
          if(ref_genome.seq[bp_tmp] != anc_genome.seq[bp_tmp]){

            //fixed snp
            if(tmrca <= outgroup_age){

              if(ref_genome.seq[bp_tmp] == 'A' || ref_genome.seq[bp_tmp] == 'C' || ref_genome.seq[bp_tmp] == 'G' || ref_genome.seq[bp_tmp] == 'T' || ref_genome.seq[bp_tmp] == '0' || ref_genome.seq[bp_tmp] == '1'){

                if(anc_genome.seq[bp_tmp] == 'A' || anc_genome.seq[bp_tmp] == 'C' || anc_genome.seq[bp_tmp] == 'G' || anc_genome.seq[bp_tmp] == 'T' || anc_genome.seq[bp_tmp] == '0' || anc_genome.seq[bp_tmp] == '1'){

                  mut_combined.info[snp_comb].branch.resize(1);
                  mut_combined.info[snp_comb].branch[0] = root;
                  mut_combined.info[snp_comb].age_begin = tmrca;
                  mut_combined.info[snp_comb].age_end   = outgroup_age;
                  assert(mut_combined.info[snp_comb].age_begin <= mut_combined.info[snp_comb].age_end);

                  mut_combined.info[snp_comb].tree = tree_count;
                  mut_combined.info[snp_comb].pos  = bp_tmp + 1;
                  mut_combined.info[snp_comb].freq.resize(1);
                  mut_combined.info[snp_comb].freq[0] = N;

                  ancestral = anc_genome.seq[bp_tmp];
                  derived   = ref_genome.seq[bp_tmp];

                  mut_combined.info[snp_comb].mutation_type = ancestral + "/" + derived;
                  if(snp_comb > 0) mut_combined.info[snp_comb-1].dist = bp_tmp + 1 - mut_combined.info[snp_comb-1].pos; 
                  mut_combined.info[snp_comb].snp_id = snp_comb;

                  snp_comb++;

                }
              }


            }

          }
        }
      }
    } 

    if(biallelic){

      //check if snp exists in mut
      if(snp_ref < L_ref){
        while((*it_mut).pos < bp+1){
          num_bases_snp_persists = ancmut.NextSNP(mtr, it_mut);
          snp_ref++;
          if(snp_ref == L_ref) break;
        }
      }

      if(tree_count < (*it_mut).tree){
        tree_count = (*it_mut).tree;
        mtr.tree.GetCoordinates(coords);
        tmrca = coords[root];
      } 

      if((*it_mut).pos == bp+1 && DAF > 0 && DAF < N){

        if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1){

          std::string ancestral, derived, vcf_reference = vcf.rec->d.allele[0], vcf_alternative = vcf.rec->d.allele[1];
          int i = 0;
          ancestral.clear();
          derived.clear();
          while((*it_mut).mutation_type[i] != '/' && i < (*it_mut).mutation_type.size()){
            ancestral.push_back((*it_mut).mutation_type[i]);
            i++;
          }
          i++;
          while(i < (*it_mut).mutation_type.size()){
            derived.push_back((*it_mut).mutation_type[i]);
            i++;
          }

          if( (ancestral == vcf_reference && derived == vcf_alternative) || (derived == vcf_reference && ancestral == vcf_alternative) ){

            if((derived == vcf_reference && ancestral == vcf_alternative)) DAF = N-DAF;

            if((*it_mut).age_end > 0){
              //if segregating, copy over
              //TODO: flag whether this is a singleton
              mut_combined.info[snp_comb] = (*it_mut);
              mut_combined.info[snp_comb].tree = tree_count;
              mut_combined.info[snp_comb].pos = bp+1;
              mut_combined.info[snp_comb].freq.clear();
              mut_combined.info[snp_comb].freq.push_back(DAF);
              if(snp_comb > 0) mut_combined.info[snp_comb-1].dist = bp+1 - mut_combined.info[snp_comb-1].pos; 
              mut_combined.info[snp_comb].snp_id = snp_comb;
              if(mut_combined.info[snp_comb].age_begin == 0.0){
                assert(DAF == 1);
              }
              snp_comb++;
            }

          }

        }

      }else{

        //SNP is not segregating in tree, so now need to check whether it is fixed for the derived allele
        //either DAF == N, vcf_alternative != anc_genome, vcf_reference == 
        //DAF == 0, and anc_reference != anc_genome
        //DAF == N, and anc_reference == anc_genome and vcf_alternative != anc_genome

        if(DAF == N || DAF == 0){

          if(tmrca <= outgroup_age && bp < ref_genome.seq.size() && bp < anc_genome.seq.size() && bp < mask.seq.size()){

            std::string vcf_reference = vcf.rec->d.allele[0], vcf_alternative = vcf.rec->d.allele[1];
            if(vcf_reference.size() == 1 && vcf_alternative.size() == 1){
              if(mask.seq[bp] == 'P'){
                if(ref_genome.seq[bp] == 'A' || ref_genome.seq[bp] == 'C' || ref_genome.seq[bp] == 'G' || ref_genome.seq[bp] == 'T' || ref_genome.seq[bp] == '0' || ref_genome.seq[bp] == '1'){

                  if(anc_genome.seq[bp] == 'A' || anc_genome.seq[bp] == 'C' || anc_genome.seq[bp] == 'G' || anc_genome.seq[bp] == 'T' || anc_genome.seq[bp] == '0' || anc_genome.seq[bp] == '1'){

                    char vcf_ref = vcf_reference[0];
                    char vcf_alt = vcf_alternative[0];

                    mut_combined.info[snp_comb].branch.resize(1);
                    mut_combined.info[snp_comb].branch[0] = root;
                    mut_combined.info[snp_comb].age_begin = tmrca;
                    mut_combined.info[snp_comb].age_end   = outgroup_age;
                    assert(mut_combined.info[snp_comb].age_begin <= mut_combined.info[snp_comb].age_end);

                    mut_combined.info[snp_comb].tree = tree_count;
                    mut_combined.info[snp_comb].pos  = bp + 1;
                    mut_combined.info[snp_comb].freq.resize(1);
                    mut_combined.info[snp_comb].freq[0] = N;


                    bool reject = false;
                    if(DAF == N){
                      if(anc_genome.seq[bp] == ref_genome.seq[bp] && vcf_ref == ref_genome.seq[bp] && vcf_alt != anc_genome.seq[bp]){
                        ancestral = anc_genome.seq[bp];
                        derived   = vcf_alt;
                      }else{
                        reject = true;
                      }
                    }else if(DAF == 0){
                      if(ref_genome.seq[bp] != anc_genome.seq[bp] && vcf_ref == ref_genome.seq[bp] && vcf_alt == anc_genome.seq[bp]){              
                        ancestral = anc_genome.seq[bp];
                        derived   = ref_genome.seq[bp];
                      }else{
                        reject = true;
                      }
                    }  

                    mut_combined.info[snp_comb].mutation_type = ancestral + "/" + derived;
                    if(snp_comb > 0) mut_combined.info[snp_comb-1].dist = bp + 1 - mut_combined.info[snp_comb-1].pos; 
                    mut_combined.info[snp_comb].snp_id = snp_comb;

                    if(reject) snp_comb--;
                    snp_comb++;

                  }
                }

              }
            }
          }

        }
      }
    }
  }

  mut_combined.info.resize(snp_comb);
  mut_combined.Dump(options["output"].as<std::string>());

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}


