#include <iostream>

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

void
mut(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, haps, sample, input, output. Optional: num_bins, coal" << std::endl;
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

  int correct = 0;
  if(options.count("correction") > 0){ 
    correct = options["correction"].as<int>();
  }

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

  }else if(1){
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

    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }else if(0){
    num_epochs = 6;
    epochs.resize(num_epochs);

    epochs[0] = 0;
    epochs[1] = 1e3/28;
    epochs[2] = 1e4/28;
    epochs[3] = 1e5/28;
    epochs[4] = 0.99e7/28;
    epochs[5] = 1e8;
  }else{
    num_epochs = 4;
    epochs.resize(num_epochs);
    epochs[0] = 0;
    epochs[1] = 1e4/28;
    epochs[2] = 1e7/28;
    epochs[3] = 1e8/28;
  }


  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)

  bool all = true;
  Data data(1,1);
  if(!all){
    haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
    data.N = ref_N.GetN();
    ref_N.CloseFile();
  }

  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);

  if(options.count("coal") > 0){
    int i = 0;
    double dummy;
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      is >> coal_rates[i][e];
      std::cerr << coal_rates[i][e] << " "; 
    }
    std::cerr << std::endl;

    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      is >> dummy >> dummy;
      for(int e = 0; e < num_epochs; e++){
        coal_rates[i][e] = coal_rates[0][e]; 
      }
    }
  }else{
    double initial_coal_rate = 1.0/20000.0;
    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }
  }

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 10;
  int num_age_bins = ((int) (log(1e8) * C))+1;
  //std::cerr << num_age_bins << std::endl;
  std::vector<double> age_bin(num_age_bins, 0.0);
  std::vector<double>::iterator it_age_bin = age_bin.begin();
  int bin = 0;
  *it_age_bin = 0.0;
  it_age_bin++;
  for(; bin < num_age_bins-1; bin++){
    *it_age_bin =  exp(bin/C)/10.0;
    it_age_bin++;
  }

  std::vector<std::vector<int>> age_shared_count(data.N), age_notshared_count(data.N);
  for(int i = 0; i < data.N; i++){
    age_shared_count[i].resize(num_age_bins*num_age_bins);
    age_notshared_count[i].resize(num_age_bins*num_age_bins);
    std::fill(age_shared_count[i].begin(), age_shared_count[i].end(), 0);
    std::fill(age_notshared_count[i].begin(), age_notshared_count[i].end(), 0);
  }

  ///////////

  std::cerr << "Parsing input files" << std::endl;

  std::uniform_real_distribution<double> dist_unif(0,1);   
  std::mt19937 rng;
  int seed = 2;
  rng.seed(seed);

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
    std::cerr << chromosomes.size() << std::endl;

  }else{
    chromosomes.resize(0);
    exit(1);
  }

  for(int chr = 0; chr < chromosomes.size(); chr++){

    std::cerr << "CHR: " << chr << std::endl;

    std::string filename_anc      = options["anc"].as<std::string>() + "_chr" + chromosomes[chr] + ".anc";
    std::string filename_mut      = options["mut"].as<std::string>() + "_chr" + chromosomes[chr] + ".mut";
    std::string filename_haps     = options["haps"].as<std::string>() + "_chr" + chromosomes[chr] + ".haps.gz";
    std::string filename_sample   = options["sample"].as<std::string>() + "_chr" + chromosomes[chr] + ".sample.gz";
    std::string file_input_haps   = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".haps.gz";
    std::string file_input_sample = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".sample.gz";

    Mutations mut;
    mut.Read(filename_mut);
    haps ref(filename_haps.c_str(), filename_sample.c_str());
    haps input(file_input_haps.c_str(), file_input_sample.c_str());

    int DAF;
    int bp_ref, bp_input = -1;
    int snp = 0, snp_input = 0, snp_ref = 0;
    double age_begin, age_end;
    std::vector<char> sequence_ref(ref.GetN()), sequence_input(input.GetN());
    //iterating over snps in reference
    for(; snp_ref < ref.GetL();){
      ref.ReadSNP(sequence_ref, bp_ref);
      snp_ref++;
      while((bp_input == -1 || bp_input < bp_ref) && snp_input < input.GetL()){
        input.ReadSNP(sequence_input, bp_input);
        snp_input++;
        if(snp_input == input.GetL()) break;
      }
      while(bp_ref > mut.info[snp].pos){
        snp++;
        if(snp == mut.info.size()) break;
      }
      if(snp == mut.info.size()) break;

      bool use_SNP = false;
      DAF = 0;
      for(std::vector<char>::iterator it_seq = sequence_ref.begin(); it_seq != sequence_ref.end(); it_seq++){
        if(*it_seq == '1') DAF++;
      }

      if(bp_ref == mut.info[snp].pos){
        //either correct == 0, or correct == 1 and DAF < N
        age_begin = mut.info[snp].age_begin;
        age_end   = mut.info[snp].age_end;
        if(mut.info[snp].flipped == 0 && mut.info[snp].branch.size() == 1 && age_begin < age_end){
          use_SNP   = true;
        }
      }

      if(use_SNP && age_end > 0){

        assert(age_begin <= age_end);

        int bin_index_age;
        double age = dist_unif(rng) * (age_end - age_begin) + age_begin;
        bin_index_age = std::max(0, (int)std::round(log(10*age)*C)+1);
        bin_index_age = bin_index_age * num_age_bins + bin_index_age;

        int bin_index1 = std::max(0, (int)std::round(log(10*age_begin)*C)+1);
        int bin_index2 = std::max(0, (int)std::round(log(10*age_end)*C)+1);

        int bin_index  = bin_index1 * num_age_bins + bin_index2;
        assert(bin_index < num_age_bins*num_age_bins);

        int count = 0;
        for(int j = 0; j < ref.GetN(); j++){
          if(sequence_ref[j] == '1'){
            for(int i = 0; i < input.GetN(); i++){
              if(bp_ref == bp_input && sequence_input[i] == '1'){
                //shared
                if(DAF == 1){
                  age_shared_count[0][bin_index]++;
                }else{
                  age_shared_count[0][bin_index_age]++;
                }
              }else{
                //not shared
                age_notshared_count[0][bin_index_age]++;
              }
            }
            count++;
          }
          if(count == DAF) break;
        }

      }
    }
    ref.CloseFile();
    input.CloseFile();

  }

  ////////////////////////////////////////  

  std::ofstream os_log(options["output"].as<std::string>() + ".log");

  int max_iter = 10000;
  int perc = -1;
  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      std::cerr << "[" << perc << "%]\r";

      std::ofstream os(options["output"].as<std::string>() + ".coal");

      for(int i = 0; i < data.N; i++){
        os << i << " ";
      }
      os << std::endl;

      for(int e = 0; e < num_epochs; e++){
        os << epochs[e] << " ";
      }
      os << std::endl;

      for(int i = 0; i < data.N; i++){
        os << "0 " << i << " ";
        for(int e = 0; e < num_epochs; e++){
          os << coal_rates[i][e] << " ";
        }
        os << std::endl;
      }

      os.close();

    }

    std::vector<coal_EM> EM;
    EM.reserve(data.N);
    for(int i = 0; i < data.N; i++){
      EM.push_back(coal_EM(epochs, coal_rates[i]));
    }

    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;

    double count = 0;
    std::vector<double> num(num_epochs,0.0), denom(num_epochs,0.0);
    for(int bin1 = 0; bin1 < num_age_bins; bin1++){
      for(int bin2 = bin1; bin2 < num_age_bins; bin2++){

        int bin = bin1*num_age_bins + bin2;
        for(int j = 0; j < data.N; j++){
          if(age_shared_count[j][bin] > 0){ 
            count = age_shared_count[j][bin];
            //if(bin1 < bin2) count /= (age_bin[bin2] - age_bin[bin1]);
            double logl = EM[j].EM_shared(age_bin[bin1], age_bin[bin2], num, denom);
            log_likelihood += count * logl;
            for(int e = 0; e < num_epochs; e++){
              assert(!std::isnan(num[e]));
              assert(!std::isnan(denom[e]));
              assert(num[e] >= 0.0);
              assert(denom[e] >= 0.0);
              coal_rates_num[j][e]   += count * num[e];
              coal_rates_denom[j][e] += count * denom[e];
            }
          }
          if(age_notshared_count[j][bin] > 0){ 
            count = age_notshared_count[j][bin];
            //if(bin1 < bin2) count /= (age_bin[bin2] - age_bin[bin1]);
            double logl = EM[j].EM_notshared(age_bin[bin1], age_bin[bin2], num, denom);
            log_likelihood += count * logl;
            for(int e = 0; e < num_epochs; e++){
              assert(!std::isnan(num[e]));
              assert(!std::isnan(denom[e]));
              assert(num[e] >= 0.0);
              assert(denom[e] >= 0.0);
              coal_rates_num[j][e]   += count * num[e];
              coal_rates_denom[j][e] += count * denom[e];
            }
          }
        }

      }	
    }

    for(int i = 0; i < data.N; i++){
      for(int e = 0; e < num_epochs; e++){
        if(coal_rates_num[i][e] == 0){
          if(e > 0){
            coal_rates[i][e] = coal_rates[i][e-1];
          }else{
            coal_rates[i][e] = 0;
          }
        }else if(coal_rates_denom[i][e] == 0){
        }else{
          coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
          if(coal_rates[i][e] < 1e-7) coal_rates[i][e] = 1e-7;
          assert(coal_rates[i][e] >= 0.0);
        }
      }
    }

    os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << std::endl;

    for(int i = 0; i < data.N; i++){
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }

  }
  os_log.close();

  std::ofstream os(options["output"].as<std::string>() + ".coal");

  for(int i = 0; i < data.N; i++){
    os << i << " ";
  }
  os << std::endl;

  for(int e = 0; e < num_epochs; e++){
    os << epochs[e] << " ";
  }
  os << std::endl;

  for(int i = 0; i < data.N; i++){
    os << "0 " << i << " ";
    for(int e = 0; e < num_epochs; e++){
      os << coal_rates[i][e] << " ";
    }
    os << std::endl;
  }

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
mut_tree(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, haps, sample, input, output. Optional: num_bins, coal" << std::endl;
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

  int correct = 0;
  if(options.count("correction") > 0){ 
    correct = options["correction"].as<int>();
  }

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

    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }

  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)

  Data data(1,1);
  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);

  if(options.count("coal") > 0){
    int i = 0;
    double dummy;
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      is >> coal_rates[i][e];
      std::cerr << coal_rates[i][e] << " "; 
    }
    std::cerr << std::endl;

    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      is >> dummy >> dummy;
      for(int e = 0; e < num_epochs; e++){
        coal_rates[i][e] = coal_rates[0][e]; 
      }
    }
  }else{
    double initial_coal_rate = 1.0/10000.0;
    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }
  }

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 100;
  int num_age_bins = ((int) (log(1e8) * C)) + 1;
  std::vector<double> age_bin(num_age_bins, 0.0);
  std::vector<double>::iterator it_age_bin = age_bin.begin();
  int bin = 0;
  *it_age_bin = 0;
  it_age_bin++;
  for(bin = 1; bin < num_age_bins; bin++){
    *it_age_bin =  exp((bin-1.0)/C)/10.0;
    it_age_bin++;
  }

  ///////////

  std::cerr << "Parsing input files" << std::endl;

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
  haps input((options["input"].as<std::string>() + ".haps.gz").c_str(), (options["input"].as<std::string>() + ".sample.gz").c_str());

  double outgroup_tmrca = 10e6/28;
  int N = ancmut.NumTips();
  int root = 2 * N - 2;
  std::vector<float> coords(root+1), coords_mut(root+1);
  std::vector<int> sorted_indices(root+1);

  std::vector<std::vector<float>> num_lin(ancmut.NumSnps());
  std::vector<std::vector<float>> DAF(ancmut.NumSnps());	
  std::vector<float> tmrca(ancmut.NumSnps(), 1e8/28.0);
  std::vector<double> age_begin(ancmut.NumSnps());
  std::vector<double> age_end(ancmut.NumSnps());
  std::vector<int> tree_index(ancmut.NumSnps());

  std::vector<std::vector<float>>::iterator it_num_lin = num_lin.begin();
  std::vector<std::vector<float>>::iterator it_DAF = DAF.begin();
  for(it_num_lin = num_lin.begin(); it_num_lin != num_lin.end(); it_num_lin++){
    (*it_num_lin).resize(num_age_bins);
  }
  for(it_DAF = DAF.begin(); it_DAF != DAF.end(); it_DAF++){
    (*it_DAF).resize(num_age_bins);
  }

  std::vector<char> sequence_input(input.GetN());
  std::vector<std::vector<int>> shared(ancmut.NumSnps());
  int N_input  = input.GetN();
  int bp_input = -1, snp_input = 0;

  //get first SNP (Only necessary when using NextSNP, for NextTree I don't need to call this)
  float num_bases_SNP_persists = ancmut.FirstSNP(mtr, it_mut);
  int tree_count = (*it_mut).tree, snp = 0;
  mtr.tree.GetCoordinates(coords);

  std::size_t m1(0);
  std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
  std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
      return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

  int lins = 0;
  double age = coords[*sorted_indices.begin()];
  int i = 0;
  std::vector<int>::iterator it_sorted_indices;
  for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    if(coords[*it_sorted_indices] > age){
      age = coords[*it_sorted_indices];
      while(age_bin[i] < age){
        num_lin[snp][i] = lins;
        i++;
        if(i == num_age_bins-1) break;
      }
    }
    if(*it_sorted_indices < N){
      lins++;
    }else{
      lins--;
    }
    assert(lins >= 1);
    if(i == num_age_bins-1) break;
  }
  for(; i < num_age_bins; i++){
    num_lin[snp][i] = 1;
  }

  //iterate through whole file
  while(num_bases_SNP_persists >= 0.0){

    if(snp == 100000) break;

    if((*it_mut).age_begin < (*it_mut).age_end && (*it_mut).age_end > 0 && (*it_mut).flipped == 0 && (*it_mut).branch.size() == 1){

      if(tree_count < (*it_mut).tree){

        tree_count = (*it_mut).tree;
        mtr.tree.GetCoordinates(coords);

        std::size_t m1(0);
        std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
        std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
            return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

        lins = 0;
        age = coords[*sorted_indices.begin()];
        i = 0;
        for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
          if(coords[*it_sorted_indices] > age){
            age = coords[*it_sorted_indices];
            while(age_bin[i] < age){
              num_lin[snp][i] = lins;
              i++;
              if(i == num_age_bins-1) break;
            }
          }
          if(*it_sorted_indices < N){
            lins++;
          }else{
            lins--;
          }
          assert(lins >= 1);
          if(i == num_age_bins-1) break;
        }
        for(; i < num_age_bins; i++){
          num_lin[snp][i] = 1;
        }

      }else{
        if(snp > 0) num_lin[snp]   = num_lin[snp-1];
      }

      while((bp_input == -1 || bp_input < (*it_mut).pos) && snp_input < input.GetL()){
        input.ReadSNP(sequence_input, bp_input);
        snp_input++;
        if(snp_input == input.GetL()) break;
      }

      shared[snp].resize(input.GetN());
      std::fill(shared[snp].begin(), shared[snp].end(), 0);
      if(bp_input == (*it_mut).pos){
        for(int j = 0; j < input.GetN(); j++){
          shared[snp][j] = (sequence_input[j] == '1');
        }
      }

      tree_index[snp] = tree_count;
      int i_end      = std::max(0.0, log((*it_mut).age_end*10)*C + 1);
      int i_begin    = std::max(0.0, log((*it_mut).age_begin*10)*C + 1);
      age_begin[snp] = age_bin[i_begin];
      age_end[snp]   = age_bin[i_end];
      //tmrca[snp]     = coords[root];

      std::fill(coords_mut.begin(), coords_mut.end(), coords[root]+1);
      std::fill(DAF[snp].begin(), DAF[snp].end(), 0);

      mtr.tree.GetCoordinates(*(*it_mut).branch.begin(), coords_mut);
      std::size_t m2(0);
      std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m2++; });
      std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
          return std::tie(coords_mut[i1],i1) < std::tie(coords_mut[i2],i2); } );

      lins = 0;
      age = coords_mut[*sorted_indices.begin()];
      i = 0;
      for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
        if(coords_mut[*it_sorted_indices] == coords[root] + 1) break;
        if(coords_mut[*it_sorted_indices] > age){
          age = coords_mut[*it_sorted_indices];
          while(age_bin[i] < age){
            DAF[snp][i] = lins;
            assert(DAF[snp][i] <= num_lin[snp][i]);
            i++;
            if(i == num_age_bins-1 || i == i_begin) break;
          }
        }
        if(*it_sorted_indices < N){
          lins++;
        }else{
          lins--;
        }
        assert(lins >= 1);
        if(i == num_age_bins-1 || i == i_begin) break;
      }
      while(age_bin[i] < age_end[snp]){
        DAF[snp][i] = 1;
        assert(DAF[snp][i] <= num_lin[snp][i]);
        i++;
        if(i == num_age_bins-1) break;
      }
      if(i_begin != i_end) {
        assert(DAF[snp][i_begin] == 1);
        if(i_begin > 0){
          assert(DAF[snp][i_begin-1] > 1);
        }
      }

      snp++;
    }

    //it_mut points to a SNP, mtr stores the marginal tree corresponding to that SNP.
    num_bases_SNP_persists = ancmut.NextSNP(mtr, it_mut);

  }
  tree_index.resize(snp);
  age_begin.resize(snp);
  age_end.resize(snp);
  DAF.resize(snp);
  num_lin.resize(snp);
  shared.resize(snp);

  ////////////////////////////////////////  

  std::ofstream os_log(options["output"].as<std::string>() + ".log");

  int max_iter = 200;
  int perc = -1;
  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    double start_time = time(NULL);
    clock_t begin = clock();

    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      //std::cerr << "[" << perc << "%]\r";

      if(0){
        std::ofstream os(options["output"].as<std::string>() + ".coal");

        for(int i = 0; i < data.N; i++){
          os << i << " ";
        }
        os << std::endl;

        for(int e = 0; e < num_epochs; e++){
          os << epochs[e] << " ";
        }
        os << std::endl;

        for(int i = 0; i < data.N; i++){
          os << "0 " << i << " ";
          for(int e = 0; e < num_epochs; e++){
            os << coal_rates[i][e] << " ";
          }
          os << std::endl;
        }

        os.close();
      }

    }

    coal_EM_tree EM(C, epochs, coal_rates[0]);

    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;

    int count = 0;
    std::vector<double> num(num_epochs, 0.0), denom(num_epochs,0.0);
    snp = 0;
    EM.UpdateTree(num_lin[snp]);
    int tree = tree_index[snp];
    for(; snp < age_begin.size(); snp++){

      if(tree < tree_index[snp]){
        EM.UpdateTree(num_lin[snp]);
        tree = tree_index[snp];
      }

      //TODO: optimise by precalculating values for each tree, precalculating values for coalescence rates
      for(int j = 0; j < N_input; j++){
        if(shared[snp][j] == 1){
          //std::cerr << snp << " shared" << std::endl;
          log_likelihood += EM.EM_shared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp], num, denom);
        }else{
          //std::cerr << snp << " not shared" << std::endl;
          log_likelihood += EM.EM_notshared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp], num, denom);
        }

        for(int e = 0; e < num_epochs; e++){
          coal_rates_num[0][e]   += num[e];
          coal_rates_denom[0][e] += denom[e];
        }
      }

    }

    for(int i = 0; i < data.N; i++){
      for(int e = 0; e < num_epochs; e++){
        //std::cerr << e << " " << coal_rates[i][e] << " " << coal_rates_num[i][e] << " " << coal_rates_denom[i][e] << std::endl;
        if(coal_rates_num[i][e] == 0){
          if(e > 0){
            coal_rates[i][e] = coal_rates[i][e-1];
          }else{
            coal_rates[i][e] = 0;
          }
        }else if(coal_rates_denom[i][e] == 0){
        }else{
          coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
          if(coal_rates[i][e] < 1e-7) coal_rates[i][e] = 1e-7;
          assert(coal_rates[i][e] >= 0.0);
        }
      }
    }

    os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << std::endl;

    for(int i = 0; i < data.N; i++){
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }

    clock_t end = clock();
    double end_time = time(NULL);
    double elapsed_secs1 = double(end - begin) / CLOCKS_PER_SEC;
    std::cerr << elapsed_secs1 << std::endl;

  }
  os_log.close();

  std::ofstream os(options["output"].as<std::string>() + ".coal");

  for(int i = 0; i < data.N; i++){
    os << i << " ";
  }
  os << std::endl;

  for(int e = 0; e < num_epochs; e++){
    os << epochs[e] << " ";
  }
  os << std::endl;

  for(int i = 0; i < data.N; i++){
    os << "0 " << i << " ";
    for(int e = 0; e < num_epochs; e++){
      os << coal_rates[i][e] << " ";
    }
    os << std::endl;
  }

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
ParseFiles(std::string& filename_anc, std::string& filename_mut, std::string& filename_haps, std::string& filename_sample, std::vector<std::vector<int>>& shared, std::vector<std::vector<float>>& num_lin, std::vector<std::vector<float>>& DAF, std::vector<double>& age_begin, std::vector<double>& age_end, std::vector<int>& tree_index, int& N, int& N_input, const int C, const std::vector<double>& age_bin, int& snp){

  std::uniform_real_distribution<double> dist_unif(0,1);   
  std::mt19937 rng;
  int seed = snp;
  rng.seed(seed);

  int num_age_bins = age_bin.size();

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(filename_anc, filename_mut); 
  haps input(filename_haps.c_str(), filename_sample.c_str());

  double outgroup_tmrca = 10e6/28;
  N = ancmut.NumTips();
  int root = 2 * N - 2;
  std::vector<float> coords(root+1), coords_mut(root+1);
  std::vector<int> sorted_indices(root+1);

  num_lin.resize(snp + ancmut.NumSnps());
  DAF.resize(snp + ancmut.NumSnps());
  //tmrca.resize(tmrca.size() + ancmut.NumSnps());
  age_begin.resize(snp + ancmut.NumSnps());
  age_end.resize(snp + ancmut.NumSnps());
  tree_index.resize(snp + ancmut.NumSnps());
  shared.resize(snp + ancmut.NumSnps());

  std::vector<char> sequence_input(input.GetN());
  N_input  = input.GetN();
  int bp_input = -1, snp_input = 0;

  //get first SNP (Only necessary when using NextSNP, for NextTree I don't need to call this)
  float num_bases_SNP_persists = ancmut.FirstSNP(mtr, it_mut);
  int tree_count = (*it_mut).tree;
  mtr.tree.GetCoordinates(coords);

  float num_muts_mapping = 0.0, threshold = 0;
  for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
    num_muts_mapping += (*it_node).num_events;
  }

  std::size_t m1(0);
  std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
  std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
      return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

  int lins = 0;
  double age = coords[*sorted_indices.begin()];
  int i = 0;
  std::vector<int>::iterator it_sorted_indices;
  num_lin[snp].resize(num_age_bins);
  for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    if(coords[*it_sorted_indices] > age){
      age = coords[*it_sorted_indices];
      while(age_bin[i] < age){
        num_lin[snp][i] = lins;
        i++;
        if(i == num_age_bins-1) break;
      }
    }
    if(*it_sorted_indices < N){
      lins++;
    }else{
      lins--;
    }
    assert(lins >= 1);
    if(i == num_age_bins-1) break;
  }
  for(; i < num_age_bins; i++){
    num_lin[snp][i] = 1;
  }

  bool first_snp = true;
  //iterate through whole file
  while(num_bases_SNP_persists >= 0.0){

    if(tree_count < (*it_mut).tree){
      num_muts_mapping = 0.0;
      for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
        num_muts_mapping += (*it_node).num_events;
      }
    }

    if((*it_mut).age_begin < (*it_mut).age_end && (*it_mut).age_end > 0 && (*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && num_muts_mapping > threshold){

      if(*(*it_mut).branch.begin() <= root){

        if(tree_count < (*it_mut).tree){

          tree_count = (*it_mut).tree;
          mtr.tree.GetCoordinates(coords);

          std::size_t m1(0);
          std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
          std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
              return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

          lins = 0;
          age = coords[*sorted_indices.begin()];
          i = 0;
          num_lin[snp].resize(num_age_bins);
          for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
            if(coords[*it_sorted_indices] > age){
              age = coords[*it_sorted_indices];
              while(age_bin[i] < age){
                num_lin[snp][i] = lins;
                i++;
                if(i == num_age_bins-1) break;
              }
            }
            if(*it_sorted_indices < N){
              lins++;
            }else{
              lins--;
            }
            assert(lins >= 1);
            if(i == num_age_bins-1) break;
          }
          for(; i < num_age_bins; i++){
            num_lin[snp][i] = 1;
          }

        }else{
          if(!first_snp){
            num_lin[snp]   = num_lin[snp-1];
          }
        }
        first_snp = false;
        assert(num_lin.size() > 0);
        //std::cerr << snp << " " << first_snp << " " << num_lin[snp].size() << std::endl;

        while((bp_input == -1 || bp_input < (*it_mut).pos) && snp_input < input.GetL()){
          input.ReadSNP(sequence_input, bp_input);
          snp_input++;
          if(snp_input == input.GetL()) break;
        }

        shared[snp].resize(input.GetN());
        std::fill(shared[snp].begin(), shared[snp].end(), 0);
        if(bp_input == (*it_mut).pos){
          for(int j = 0; j < input.GetN(); j++){
            shared[snp][j] = (sequence_input[j] == '1');
          }
        }

        tree_index[snp] = tree_count;
        int i_end      = std::max(0.0, log((*it_mut).age_end*10)*C + 1.0);
        int i_begin    = std::max(0.0, log((*it_mut).age_begin*10)*C + 1.0);
        if(1){	
          age_begin[snp] = age_bin[i_begin];
          age_end[snp]   = age_bin[i_end];
        }else{
          double age = dist_unif(rng) * ((*it_mut).age_end - (*it_mut).age_begin) + (*it_mut).age_begin;
          int i = std::max(0.0, log(age*10)*C + 1.0);
          age_begin[snp] = age_bin[i];
          age_end[snp] = age_bin[i];
        }
        //tmrca[snp]     = coords[root];

        std::fill(coords_mut.begin(), coords_mut.end(), coords[root]+1);
        DAF[snp].resize(num_age_bins);
        std::fill(DAF[snp].begin(), DAF[snp].end(), 0);

        assert(*(*it_mut).branch.begin() <= root);
        mtr.tree.GetCoordinates(*(*it_mut).branch.begin(), coords_mut);
        std::size_t m2(0);
        std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m2++; });
        std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
            return std::tie(coords_mut[i1],i1) < std::tie(coords_mut[i2],i2); } );

        lins = 0;
        age = coords_mut[*sorted_indices.begin()];
        i = 0;
        for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
          if(coords_mut[*it_sorted_indices] == coords[root] + 1) break;
          if(coords_mut[*it_sorted_indices] > age){
            age = coords_mut[*it_sorted_indices];
            while(age_bin[i] < age){
              DAF[snp][i] = lins;
              assert(DAF[snp][i] <= num_lin[snp][i]);
              DAF[snp][i] /= num_lin[snp][i];
              i++;
              if(i == num_age_bins-1 || i == i_begin) break;
            }
          }
          if(*it_sorted_indices < N){
            lins++;
          }else{
            lins--;
          }
          assert(lins >= 1);
          if(i == num_age_bins-1 || i == i_begin) break;
        }
        while(age_bin[i] < age_end[snp]){
          DAF[snp][i] = 1;
          assert(DAF[snp][i] <= num_lin[snp][i]);
          DAF[snp][i] /= num_lin[snp][i];
          i++;
          if(i == num_age_bins-1) break;
        }

        snp++;

      }
    }

    //it_mut points to a SNP, mtr stores the marginal tree corresponding to that SNP.
    num_bases_SNP_persists = ancmut.NextSNP(mtr, it_mut);

  }

}

void
mut_tree_fast(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("mut") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, input, output. Optional: num_bins, coal" << std::endl;
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

  int regularise = 0;
  if(options.count("regularise") > 0){ 
    regularise = options["regularise"].as<int>();
  }

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

  }else if(1){
    num_epochs = 20;
    if(options.count("num_bins") > 0){
      num_epochs = options["num_bins"].as<int>();
    }
    float years_per_gen = 28.0;
    if(options.count("years_per_gen")){
      years_per_gen = options["years_per_gen"].as<float>();
    }
    num_epochs++; 
    epochs.resize(num_epochs);

    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;
  }else if(0){
    num_epochs = 6;
    epochs.resize(num_epochs);

    epochs[0] = 0;
    epochs[1] = 1e3/28;
    epochs[2] = 1e4/28;
    epochs[3] = 1e5/28;
    epochs[4] = 0.99e7/28;
    epochs[5] = 1e8;
  }else{
    num_epochs = 4;
    epochs.resize(num_epochs);
    epochs[0] = 0;
    epochs[1] = 1e4/28;
    epochs[2] = 1e7/28;
    epochs[3] = 1e8/28;
  }

  std::cerr << num_epochs << std::endl;
  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)

  Data data(1,1);
  //std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);
  std::vector<double> coal_rates(num_epochs), coal_rates_num(num_epochs), coal_rates_denom(num_epochs);

  if(options.count("coal") > 0){
    double dummy;
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      is >> coal_rates[e];
      std::cerr << coal_rates[e] << " "; 
    }
    std::cerr << std::endl;
  }else{
    std::fill(coal_rates.begin(), coal_rates.end(), 0.0);
    std::fill(coal_rates_num.begin(), coal_rates_num.end(), 0.0);
    std::fill(coal_rates_denom.begin(), coal_rates_denom.end(), 0.0);
  }

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 500;
  int num_age_bins = ((int) (log(1e8) * C)) + 1;
  std::vector<double> age_bin(num_age_bins, 0.0);
  std::vector<double>::iterator it_age_bin = age_bin.begin();
  int bin = 0;
  *it_age_bin = 0;
  it_age_bin++;
  for(bin = 1; bin < num_age_bins; bin++){
    *it_age_bin =  exp((bin-1.0)/C)/10.0;
    it_age_bin++;
  }

  ///////////

  std::cerr << "Parsing input files" << std::endl;

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
    std::cerr << chromosomes.size() << std::endl;

  }else{
    chromosomes.resize(0);
  }

  std::vector<std::vector<int>> shared;
  std::vector<std::vector<float>> num_lin;
  std::vector<std::vector<float>> DAF;	
  std::vector<double> age_begin;
  std::vector<double> age_end;
  std::vector<int> tree_index;
  int snp = 0;
  int N, N_input;

  if(chromosomes.size() > 0){
    for(int chr = 0; chr < chromosomes.size(); chr++){
      std::cerr << "CHR: " << chromosomes[chr] << std::endl;
      std::string filename_anc    = options["anc"].as<std::string>() + "_chr" + chromosomes[chr] + ".anc";
      std::string filename_mut    = options["mut"].as<std::string>() + "_chr" + chromosomes[chr] + ".mut";
      std::string filename_haps   = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".haps.gz";
      std::string filename_sample = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".sample.gz";
      ParseFiles(filename_anc, filename_mut, filename_haps, filename_sample, shared, num_lin, DAF, age_begin, age_end, tree_index, N, N_input, C, age_bin, snp);
    }
  }else{
    std::string filename_anc    = options["anc"].as<std::string>();
    std::string filename_mut    = options["mut"].as<std::string>();
    std::string filename_haps   = options["input"].as<std::string>() + ".haps.gz";
    std::string filename_sample = options["input"].as<std::string>() + ".sample.gz";
    ParseFiles(filename_anc, filename_mut, filename_haps, filename_sample, shared, num_lin, DAF, age_begin, age_end, tree_index, N, N_input, C, age_bin, snp);
  }
  tree_index.resize(snp);
  age_begin.resize(snp);
  age_end.resize(snp);
  DAF.resize(snp);
  num_lin.resize(snp);
  shared.resize(snp);

  std::cerr << "Using " << snp << " SNPs" << std::endl << std::endl;

  ////////////////////////////////////////  

  std::ofstream os_log(options["output"].as<std::string>() + ".log");

  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);

  coal_EM_tree_fast EM(C, N, epochs);
  if(options.count("coal") == 0){
    std::vector<double> epochs_tmp(2), coal_rates_tmp(2);
    epochs_tmp[0] = 0;
    epochs_tmp[1] = 1e8/28;
    EM.UpdateEpochs(epochs_tmp);

    double coal_rate_min = 1e-6, coal_rate;
    double logl_max = log(0.0);
    int arg_max = 1;
    if(1){
      for(int f = 1; f <= 10; f++){

        coal_rate = coal_rate_min * exp(log(10) * 0.3*(f-1));
        for(int e = 0; e < epochs_tmp.size(); e++){
          coal_rates_tmp[e] = coal_rate;
        }

        //initialise EM
        EM.UpdateCoal(coal_rates_tmp);
        snp = 0;
        log_likelihood = 0.0;
        EM.UpdateTree(num_lin[snp]);
        int tree = tree_index[snp];
        for(; snp < age_begin.size(); snp++){

          if(tree != tree_index[snp]){
            EM.UpdateTree(num_lin[snp]);
            tree = tree_index[snp];
          }

          bool is_shared = false, is_notshared = false;	
          double logl_shared, logl_notshared;
          for(int j = 0; j < N_input; j++){
            if(shared[snp][j] == 1){
              if(is_shared == false){
                logl_shared = EM.Logl_shared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp]);
                is_shared = true;
              }
              log_likelihood += logl_shared;
            }else{
              if(is_notshared == false){
                logl_notshared = EM.Logl_notshared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp]);
                is_notshared = true;
              }
              log_likelihood += logl_notshared;
            }
          }

        }

        std::cerr << f << " " << 0.5/coal_rate << " " << log_likelihood << std::endl;
        if(logl_max < log_likelihood){
          logl_max = log_likelihood;
          arg_max  = f;
        }

      }
    }
    coal_rate = coal_rate_min * exp(log(10) * 0.3*(arg_max-1));
    std::fill(coal_rates.begin(), coal_rates.end(), coal_rate);

    std::cerr << arg_max << " " << 0.5/coal_rate << std::endl;
    std::cerr << std::endl;

  }

  double gamma = 1;
  std::vector<double> gradient(num_epochs, 0.0), gradient_prev(num_epochs, 0.0), coal_rates_prev(num_epochs, 0.0), coal_rates_nesterov(num_epochs, 0.0);
  coal_rates_nesterov = coal_rates;

  int max_iter = 1000;
  if(regularise == 2){
    max_iter = 5000;
  }
  int perc = -1;

  double start_time = time(NULL);
  clock_t begin = clock();

  EM.UpdateEpochs(epochs);
  log_likelihood = log(0.0);
  prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      std::cerr << "[" << perc << "%]\r";

      if(1){
        std::ofstream os(options["output"].as<std::string>() + ".coal");

        os << "0\n";

        for(int e = 0; e < num_epochs; e++){
          os << epochs[e] << " ";
        }
        os << std::endl;

        os << "0 0" << " ";
        for(int e = 0; e < num_epochs; e++){
          os << coal_rates[e] << " ";
        }
        os << std::endl;

        os.close();
      }

    }

    EM.UpdateCoal(coal_rates_nesterov);

    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;

    int count = 0;
    std::vector<double> num(num_epochs, 0.0), denom(num_epochs,0.0), num2(num_epochs, 0.0), denom2(num_epochs, 0.0);
    snp = 0;
    EM.UpdateTree(num_lin[snp]);

    int tree = tree_index[snp];
    for(; snp < age_begin.size(); snp++){

      if(tree != tree_index[snp]){
        EM.UpdateTree(num_lin[snp]);
        tree = tree_index[snp];
      }

      bool is_shared = false, is_notshared = false;
      double logl_shared, logl_notshared;
      for(int j = 0; j < N_input; j++){
        if(shared[snp][j] == 1){
          //std::cerr << snp << " shared" << std::endl;
          if(is_shared == false){
            logl_shared = EM.EM_shared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp], num, denom);
            is_shared = true;
          }
          log_likelihood += logl_shared;
          for(int e = 0; e < num_epochs; e++){
            coal_rates_num[e]   += num[e];
            coal_rates_denom[e] += denom[e];
          }
        }else{
          //std::cerr << snp << " not shared" << std::endl;
          if(is_notshared == false){
            logl_notshared = EM.EM_notshared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp], num2, denom2);
            is_notshared = true;
          }
          log_likelihood += logl_notshared;
          for(int e = 0; e < num_epochs; e++){
            coal_rates_num[e]   += num2[e];
            coal_rates_denom[e] += denom2[e];
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
      gamma = 1e-7/age_begin.size();
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
          //gamma         = coal_rates[e]/coal_rates_denom[e];
          coal_rates[e]  += 0.9*gradient_prev[e] + gamma * gradient[e];
          coal_rates_nesterov[e] = coal_rates[e] + 0.9*(0.9*gradient_prev[e] + gamma * gradient[e]);
          //coal_rates[e] += gamma * gradient[e]; 
          //coal_rates_nesterov[e] = coal_rates[e];
        }
        if(coal_rates[e] < 1e-7){
          coal_rates[e] = 1e-7;
        }
        if(coal_rates_nesterov[e] < 1e-7){
          coal_rates_nesterov[e] = 1e-7;
        }
        assert(coal_rates[e] >= 0.0);
        assert(coal_rates_nesterov[e] >= 0.0);

      }

      if(!is_EM){
        //gradient_prev[e]   = gradient[e];
        gradient_prev[e] = 0.9*gradient_prev[e] + gamma * gradient[e];
      }

    }

    os_log << gamma << " " << log_likelihood << " " << log_likelihood - prev_log_likelihood << std::endl;

    std::fill(coal_rates_num.begin(), coal_rates_num.end(), 0.0);
    std::fill(coal_rates_denom.begin(), coal_rates_denom.end(), 0.0);

  }
  os_log.close();

  clock_t end = clock();
  double end_time = time(NULL);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cerr << "Iteration time: " << elapsed_secs << std::endl;

  std::ofstream os(options["output"].as<std::string>() + ".coal");

  os << "0\n";
  for(int e = 0; e < num_epochs; e++){
    os << epochs[e] << " ";
  }
  os << std::endl;

  os << "0 0" << " ";
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
mut_perturb(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("coal") || !options.count("anc") || !options.count("mut") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: coal, anc, mut, input, output." << std::endl;
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

  int correct = 0;
  if(options.count("correction") > 0){ 
    correct = options["correction"].as<int>();
  }

  ////////////////////////////////////////

  //decide on epochs
  int num_epochs = 0; 
  std::vector<double> epochs;

  std::ifstream is;
  std::string line;
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

  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)

  Data data(1,1);
  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);

  int i = 0;
  double dummy;
  coal_rates[i].resize(num_epochs);
  coal_rates_num[i].resize(num_epochs);
  coal_rates_denom[i].resize(num_epochs);
  is >> dummy >> dummy;
  for(int e = 0; e < num_epochs; e++){
    is >> coal_rates[i][e];
    std::cerr << coal_rates[i][e] << " "; 
  }
  std::cerr << std::endl;

  for(int i = 0; i < data.N; i++){
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      coal_rates[i][e] = coal_rates[0][e]; 
    }
  }

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 5;
  int num_age_bins = ((int) (log(1e8) * C)) + 1;
  std::vector<double> age_bin(num_age_bins, 0.0);
  std::vector<double>::iterator it_age_bin = age_bin.begin();
  int bin = 0;
  *it_age_bin = 0;
  it_age_bin++;
  for(bin = 1; bin < num_age_bins; bin++){
    *it_age_bin =  exp((bin-1.0)/C)/10.0;
    it_age_bin++;
  }

  ///////////

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
    std::cerr << chromosomes.size() << std::endl;

  }else{
    chromosomes.resize(1);
    chromosomes[0] = "1";
  }

  std::cerr << "Parsing input files" << std::endl;

  std::vector<std::vector<int>> shared;
  std::vector<std::vector<float>> num_lin;
  std::vector<std::vector<float>> DAF;	
  std::vector<double> age_begin;
  std::vector<double> age_end;
  std::vector<int> tree_index;
  int snp = 0;
  int N, N_input;

  for(int chr = 0; chr < chromosomes.size(); chr++){

    std::cerr << "CHR: " << chromosomes[chr] << std::endl;

    //std::string filename_anc    = "chr_" + chromosomes[chr] + "/" + options["anc"].as<std::string>(); 
    //std::string filename_mut    = "chr_" + chromosomes[chr] + "/" + options["mut"].as<std::string>();
    //std::string filename_haps   = "chr_" + chromosomes[chr] + "/" + options["input"].as<std::string>() + ".haps.gz";
    //std::string filename_sample = "chr_" + chromosomes[chr] + "/" + options["input"].as<std::string>() + ".sample.gz";
    std::string filename_anc    = options["anc"].as<std::string>() + "_chr" + chromosomes[chr] + ".anc";
    std::string filename_mut    = options["mut"].as<std::string>() + "_chr" + chromosomes[chr] + ".mut";
    std::string filename_haps   = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".haps.gz";
    std::string filename_sample = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".sample.gz";
    ParseFiles(filename_anc, filename_mut, filename_haps, filename_sample, shared, num_lin, DAF, age_begin, age_end, tree_index, N, N_input, C, age_bin, snp);

  }
  tree_index.resize(snp);
  age_begin.resize(snp);
  age_end.resize(snp);
  DAF.resize(snp);
  num_lin.resize(snp);
  shared.resize(snp);

  std::cerr << "Using " << snp << " SNPs" << std::endl << std::endl;

  ////////////////////////////////////////  

  std::ofstream os(options["output"].as<std::string>() + ".log");

  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);

  std::vector<double> coal_rates_tmp(num_epochs);
  coal_EM_tree_fast EM(C, N, epochs);

  int epoch_index = 0;
  double coal_rate_min = 1e-6, coal_rate1, coal_rate2;
  double logl_max = log(0.0);
  int arg_max1 = 1, arg_max2 = 1;
  for(int f = 0; f <= 20; f++){

    coal_rate1 = 0.5/(5000 + (50000 - 5000) * (f/20.0));
    coal_rates_tmp = coal_rates[0];
    coal_rates_tmp[epoch_index] = coal_rate1;

    //initialise EM
    EM.UpdateCoal(coal_rates_tmp);
    snp = 0;
    log_likelihood = 0.0;
    EM.UpdateTree(num_lin[snp]);
    int tree = tree_index[snp];
    for(; snp < age_begin.size(); snp++){

      if(tree != tree_index[snp]){
        EM.UpdateTree(num_lin[snp]);
        tree = tree_index[snp];
      }

      bool is_shared = false, is_notshared = false;	
      double logl_shared, logl_notshared;
      for(int j = 0; j < N_input; j++){
        if(shared[snp][j] == 1){
          if(is_shared == false){
            logl_shared = EM.Logl_shared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp]);
            is_shared = true;
          }
          log_likelihood += logl_shared;
        }else{
          if(is_notshared == false){
            logl_notshared = EM.Logl_notshared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp]);
            is_notshared = true;
          }
          log_likelihood += logl_notshared;
        }
      }

    }

    os << f << " " << 0.5/coal_rate1 << " " << std::setprecision(9) << log_likelihood << std::endl;
    if(logl_max < log_likelihood){
      logl_max = log_likelihood;
      arg_max1  = f;
    }

  }
  //coal_rate1 = coal_rate_min * exp(log(10) * 0.05*(arg_max1-1));
  //coal_rate2 = coal_rate_min * exp(log(10) * 0.05*(arg_max2-1));

  //std::cerr << arg_max1 << " " << 0.5/coal_rate1 << std::endl;
  //std::cerr << arg_max2 << " " << 0.5/coal_rate2 << std::endl;
  //std::cerr << std::endl;

  os.close();

  /*
     std::ofstream os(options["output"].as<std::string>() + ".coal");

     for(int i = 0; i < data.N; i++){
     os << i << " ";
     }
     os << std::endl;

     for(int e = 0; e < num_epochs; e++){
     os << epochs[e] << " ";
     }
     os << std::endl;

     for(int i = 0; i < data.N; i++){
     os << "0 " << i << " ";
     for(int e = 0; e < num_epochs; e++){
     os << coal_rates[i][e] << " ";
     }
     os << std::endl;
     }

     os.close();
     */

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
mut_logl(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("mut") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, input, output. Optional: num_bins, coal" << std::endl;
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

  int correct = 0;
  if(options.count("correction") > 0){ 
    correct = options["correction"].as<int>();
  }

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

  }else if(1){
    num_epochs = 20;
    if(options.count("num_bins") > 0){
      num_epochs = options["num_bins"].as<int>();
    }
    float years_per_gen = 28.0;
    if(options.count("years_per_gen")){
      years_per_gen = options["years_per_gen"].as<float>();
    }
    num_epochs++; 
    epochs.resize(num_epochs);

    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }else{

    num_epochs = 3;
    float years_per_gen = 28.0;
    if(options.count("years_per_gen")){
      years_per_gen = options["years_per_gen"].as<float>();
    }
    num_epochs++; 
    epochs.resize(num_epochs);

    epochs[0] = 0.0;
    epochs[1] = 5e4/28;
    epochs[2] = 1e7/28;
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }

  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)

  Data data(1,1);
  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);

  if(options.count("coal") > 0){
    int i = 0;
    double dummy;
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      is >> coal_rates[i][e];
      std::cerr << coal_rates[i][e] << " "; 
    }
    std::cerr << std::endl;

    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      is >> dummy >> dummy;
      for(int e = 0; e < num_epochs; e++){
        coal_rates[i][e] = coal_rates[0][e]; 
      }
    }
  }else{
    double initial_coal_rate = 1.0/10000.0;
    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }
  }

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 5;
  int num_age_bins = ((int) (log(1e8) * C)) + 1;
  std::vector<double> age_bin(num_age_bins, 0.0);
  std::vector<double>::iterator it_age_bin = age_bin.begin();
  int bin = 0;
  *it_age_bin = 0;
  it_age_bin++;
  for(bin = 1; bin < num_age_bins; bin++){
    *it_age_bin =  exp((bin-1.0)/C)/10.0;
    it_age_bin++;
  }

  ///////////

  std::cerr << "Parsing input files" << std::endl;

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
  haps input((options["input"].as<std::string>() + ".haps.gz").c_str(), (options["input"].as<std::string>() + ".sample.gz").c_str());

  double outgroup_tmrca = 10e6/28;
  int N = ancmut.NumTips();
  int root = 2 * N - 2;
  std::vector<float> coords(root+1), coords_mut(root+1);
  std::vector<int> sorted_indices(root+1);

  std::vector<std::vector<float>> num_lin(ancmut.NumSnps());
  std::vector<std::vector<float>> DAF(ancmut.NumSnps());	
  std::vector<float> tmrca(ancmut.NumSnps(), 1e8/28.0);
  std::vector<double> age_begin(ancmut.NumSnps());
  std::vector<double> age_end(ancmut.NumSnps());
  std::vector<int> tree_index(ancmut.NumSnps());

  std::vector<std::vector<float>>::iterator it_num_lin = num_lin.begin();
  std::vector<std::vector<float>>::iterator it_DAF = DAF.begin();
  for(it_num_lin = num_lin.begin(); it_num_lin != num_lin.end(); it_num_lin++){
    (*it_num_lin).resize(num_age_bins);
  }
  for(it_DAF = DAF.begin(); it_DAF != DAF.end(); it_DAF++){
    (*it_DAF).resize(num_age_bins);
  }

  std::vector<char> sequence_input(input.GetN());
  std::vector<std::vector<int>> shared(ancmut.NumSnps());
  int N_input  = input.GetN();
  int bp_input = -1, snp_input = 0;

  //get first SNP (Only necessary when using NextSNP, for NextTree I don't need to call this)
  float num_bases_SNP_persists = ancmut.FirstSNP(mtr, it_mut);
  int tree_count = (*it_mut).tree, snp = 0;
  mtr.tree.GetCoordinates(coords);

  std::size_t m1(0);
  std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
  std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
      return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

  int lins = 0;
  double age = coords[*sorted_indices.begin()];
  int i = 0;
  std::vector<int>::iterator it_sorted_indices;
  for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    if(coords[*it_sorted_indices] > age){
      age = coords[*it_sorted_indices];
      while(age_bin[i] < age){
        num_lin[snp][i] = lins;
        i++;
        if(i == num_age_bins-1) break;
      }
    }
    if(*it_sorted_indices < N){
      lins++;
    }else{
      lins--;
    }
    assert(lins >= 1);
    if(i == num_age_bins-1) break;
  }
  for(; i < num_age_bins; i++){
    num_lin[snp][i] = 1;
  }

  //iterate through whole file
  while(num_bases_SNP_persists >= 0.0){

    if((*it_mut).age_begin < (*it_mut).age_end && (*it_mut).age_end > 0 && (*it_mut).flipped == 0 && (*it_mut).branch.size() == 1){

      if(tree_count < (*it_mut).tree){

        tree_count = (*it_mut).tree;
        mtr.tree.GetCoordinates(coords);

        std::size_t m1(0);
        std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
        std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
            return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

        lins = 0;
        age = coords[*sorted_indices.begin()];
        i = 0;
        for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
          if(coords[*it_sorted_indices] > age){
            age = coords[*it_sorted_indices];
            while(age_bin[i] < age){
              num_lin[snp][i] = lins;
              i++;
              if(i == num_age_bins-1) break;
            }
          }
          if(*it_sorted_indices < N){
            lins++;
          }else{
            lins--;
          }
          assert(lins >= 1);
          if(i == num_age_bins-1) break;
        }
        for(; i < num_age_bins; i++){
          num_lin[snp][i] = 1;
        }

      }else{
        if(snp > 0) num_lin[snp]   = num_lin[snp-1];
      }

      while((bp_input == -1 || bp_input < (*it_mut).pos) && snp_input < input.GetL()){
        input.ReadSNP(sequence_input, bp_input);
        snp_input++;
        if(snp_input == input.GetL()) break;
      }

      shared[snp].resize(input.GetN());
      std::fill(shared[snp].begin(), shared[snp].end(), 0);
      if(bp_input == (*it_mut).pos){
        for(int j = 0; j < input.GetN(); j++){
          shared[snp][j] = (sequence_input[j] == '1');
        }
      }

      tree_index[snp] = tree_count;
      int i_end      = std::max(0.0, log((*it_mut).age_end*10)*C + 1);
      int i_begin    = std::max(0.0, log((*it_mut).age_begin*10)*C + 1);
      age_begin[snp] = age_bin[i_begin];
      age_end[snp]   = age_bin[i_end];
      //tmrca[snp]     = coords[root];

      std::fill(coords_mut.begin(), coords_mut.end(), coords[root]+1);
      std::fill(DAF[snp].begin(), DAF[snp].end(), 0);

      mtr.tree.GetCoordinates(*(*it_mut).branch.begin(), coords_mut);
      std::size_t m2(0);
      std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m2++; });
      std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
          return std::tie(coords_mut[i1],i1) < std::tie(coords_mut[i2],i2); } );

      lins = 0;
      age = coords_mut[*sorted_indices.begin()];
      i = 0;
      for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
        if(coords_mut[*it_sorted_indices] == coords[root] + 1) break;
        if(coords_mut[*it_sorted_indices] > age){
          age = coords_mut[*it_sorted_indices];
          while(age_bin[i] < age){
            DAF[snp][i] = lins;
            assert(DAF[snp][i] <= num_lin[snp][i]);
            DAF[snp][i] /= num_lin[snp][i];
            i++;
            if(i == num_age_bins-1 || i == i_begin) break;
          }
        }
        if(*it_sorted_indices < N){
          lins++;
        }else{
          lins--;
        }
        assert(lins >= 1);
        if(i == num_age_bins-1 || i == i_begin) break;
      }
      while(age_bin[i] < age_end[snp]){
        DAF[snp][i] = 1;
        assert(DAF[snp][i] <= num_lin[snp][i]);
        DAF[snp][i] /= num_lin[snp][i];
        i++;
        if(i == num_age_bins-1) break;
      }

      snp++;
    }

    //it_mut points to a SNP, mtr stores the marginal tree corresponding to that SNP.
    num_bases_SNP_persists = ancmut.NextSNP(mtr, it_mut);

  }
  tree_index.resize(snp);
  age_begin.resize(snp);
  age_end.resize(snp);
  DAF.resize(snp);
  num_lin.resize(snp);
  shared.resize(snp);

  ////////////////////////////////////////  

  std::ofstream os(options["output"].as<std::string>() + ".log");

  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);

  num_epochs = 3;
  std::vector<double> epochs_tmp(num_epochs), coal_rates_tmp(num_epochs);
  epochs_tmp[0] = 0;
  epochs_tmp[1] = 1e4/28;
  epochs_tmp[2] = 1e8/28;
  coal_EM_tree_fast EM(C, N, epochs_tmp);

  double coal_rate_min = 1e-6, coal_rate1, coal_rate2;
  double logl_max = log(0.0);
  int arg_max1 = 1, arg_max2 = 1;
  for(int f = 0; f <= 20; f++){

    for(int g = 0; g <= 20; g++){

      //coal_rate1 = coal_rate_min * exp(log(10) * 0.3*(f-1));
      //coal_rate2 = coal_rate_min * exp(log(10) * 0.3*(g-1));
      coal_rate1 = 0.5/(5000 + (40000 - 5000) * (f/20.0));
      coal_rate2 = 0.5/(5000 + (40000 - 5000) * (g/20.0));
      coal_rates_tmp[0] = coal_rate1;
      coal_rates_tmp[1] = coal_rate2;
      coal_rates_tmp[2] = coal_rate2;

      //initialise EM
      EM.UpdateCoal(coal_rates_tmp);
      snp = 0;
      log_likelihood = 0.0;
      EM.UpdateTree(num_lin[snp]);
      int tree = tree_index[snp];
      for(; snp < age_begin.size(); snp++){

        if(tree < tree_index[snp]){
          EM.UpdateTree(num_lin[snp]);
          tree = tree_index[snp];
        }

        bool is_shared = false, is_notshared = false;	
        double logl_shared, logl_notshared;
        for(int j = 0; j < N_input; j++){
          if(shared[snp][j] == 1){
            if(is_shared == false){
              logl_shared = EM.Logl_shared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp]);
              is_shared = true;
            }
            log_likelihood += logl_shared;
          }else{
            if(is_notshared == false){
              logl_notshared = EM.Logl_notshared(age_begin[snp], age_end[snp], num_lin[snp], DAF[snp]);
              is_notshared = true;
            }
            log_likelihood += logl_notshared;
          }
        }

      }

      os << f << " " << g << " " << 0.5/coal_rate1 << " " << 0.5/coal_rate2 << " " << log_likelihood << std::endl;
      if(logl_max < log_likelihood){
        logl_max = log_likelihood;
        arg_max1  = f;
        arg_max2  = g;
      }

    }

  }
  //coal_rate1 = coal_rate_min * exp(log(10) * 0.05*(arg_max1-1));
  //coal_rate2 = coal_rate_min * exp(log(10) * 0.05*(arg_max2-1));

  //std::cerr << arg_max1 << " " << 0.5/coal_rate1 << std::endl;
  //std::cerr << arg_max2 << " " << 0.5/coal_rate2 << std::endl;
  //std::cerr << std::endl;

  os.close();

  /*
     std::ofstream os(options["output"].as<std::string>() + ".coal");

     for(int i = 0; i < data.N; i++){
     os << i << " ";
     }
     os << std::endl;

     for(int e = 0; e < num_epochs; e++){
     os << epochs[e] << " ";
     }
     os << std::endl;

     for(int i = 0; i < data.N; i++){
     os << "0 " << i << " ";
     for(int e = 0; e < num_epochs; e++){
     os << coal_rates[i][e] << " ";
     }
     os << std::endl;
     }

     os.close();
     */

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

////////////////////////////////////////

void
mut_fast_simplified(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, haps, sample, input, output. Optional: num_bins, coal" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Calculating coalescence rates for (ancient) sample.." << std::endl;

  /////////////////////////

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());

  //assign an actual date to each SNP
  std::uniform_real_distribution<double> dist_unif(0,1);   
  std::mt19937 rng;
  std::vector<double> age(mut.info.size(), 0.0);
  std::vector<double>::iterator it_age = age.begin();
  int seed = 2;
  rng.seed(seed);
  for(Muts::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++ ){
    *it_age = dist_unif(rng) * ((*it_mut).age_end - (*it_mut).age_begin) + (*it_mut).age_begin;
    it_age++;
  }

  /////////////////////////////////
  //get TMRCA at each SNP
  std::vector<float> tmrca(mut.info.size(), 1e8/28.0);
  std::vector<float>::iterator it_tmrca = tmrca.begin();

  bool correct = 0;
  if(options.count("correction") > 0){ 
    correct = options["correction"].as<int>();
  }

  if(correct){
    MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
    Muts::iterator it_mut; //iterator for mut file
    AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());

    float num_bases_tree_persists = 0.0;
    int root = 2 * ancmut.NumTips() - 2;
    std::vector<float> coords(2*ancmut.NumTips() - 1);
    int tree_count = 0;
    while(num_bases_tree_persists >= 0.0){
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      if(num_bases_tree_persists >= 0.0 && (*it_mut).tree == tree_count){
        mtr.tree.GetCoordinates(coords);
        int tree_index = (*it_mut).tree;
        while((*it_mut).tree == tree_index){
          *it_tmrca = coords[root];
          it_tmrca++;
          it_mut++;
          if(it_mut == ancmut.mut_end()) break;
        }
      }
      tree_count++;
    }
    ancmut.CloseFiles();
    std::ofstream os_tmrca(options["input"].as<std::string>() + ".tmrca");
    for(it_tmrca = tmrca.begin(); it_tmrca != tmrca.end(); it_tmrca++){
      os_tmrca << *it_tmrca << "\n";
    }
    os_tmrca.close();
  }

  if(correct){
    std::ifstream is_tmrca(options["input"].as<std::string>() + ".tmrca");
    for(it_tmrca = tmrca.begin(); it_tmrca != tmrca.end(); it_tmrca++){
      is_tmrca >> *it_tmrca;
    }
    is_tmrca.close();
  }

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
      assert(epochs[e] > epochs[e-1]);
    }

  }else{
    num_epochs = 20;
    if(options.count("num_bins") > 0){
      num_epochs = options["num_bins"].as<int>();
    }
    float years_per_gen = 28.0;
    if(options.count("years_per_gen")){
      years_per_gen = options["years_per_gen"].as<float>();
    }
    num_epochs++; 
    epochs.resize(num_epochs);

    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }


  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)
  //For now, assume both are haps/sample and span same positions

  haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  Data data(ref_N.GetN(), 1);
  ref_N.CloseFile();
  //Data data(1,1);
  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);

  if(options.count("coal") > 0){
    int i = 0;
    double dummy;
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      is >> coal_rates[i][e];
    }

    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      is >> dummy >> dummy;
      for(int e = 0; e < num_epochs; e++){
        coal_rates[i][e] = coal_rates[0][e]; 
      }
    }
    is.close();
  }else{
    double initial_coal_rate = 1.0/10000.0;
    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }
  }

  ////////////////////////////////////////  

  std::ofstream os_log(options["output"].as<std::string>() + ".log");

  //for every input & ref haplotype pair, want to know where they share and not share
  //for ref haplotype, record age of SNPs that share & age of SNPs that don't share
  //const*N_ref*N_input
  //could even bin this by discretised age if wanted

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 1e2;
  int num_age_bins = ((int) (log(1e8) * C));
  //std::cerr << num_age_bins << std::endl;
  std::vector<double> age_bin(num_age_bins, 0.0);
  int bin = 0;
  for(std::vector<double>::iterator it_age_bin = age_bin.begin(); it_age_bin != age_bin.end(); it_age_bin++){
    *it_age_bin = exp(bin/C)/10.0;
    bin++;
  }

  std::vector<std::vector<int>> age_shared_count(data.N), age_notshared_count(data.N);
  for(int i = 0; i < data.N; i++){
    age_shared_count[i].resize(num_age_bins);
    age_notshared_count[i].resize(num_age_bins);
    std::fill(age_shared_count[i].begin(), age_shared_count[i].end(), 0);
    std::fill(age_notshared_count[i].begin(), age_notshared_count[i].end(), 0);
  }

  ///////////

  std::cerr << "Parsing input files" << std::endl;

  haps ref(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  haps input((options["input"].as<std::string>() + ".haps.gz").c_str(), (options["input"].as<std::string>() + ".sample.gz").c_str());

  int bp_ref, bp_input = -1;
  int snp = 0, snp_input = 0, snp_ref = 0;
  std::vector<char> sequence_ref(ref.GetN()), sequence_input(input.GetN());
  //iterating over snps in reference
  int num_used_snps = 0;
  for(; snp_ref < ref.GetL();){
    ref.ReadSNP(sequence_ref, bp_ref);
    //std::cerr << bp_ref << std::endl;
    snp_ref++;
    while((bp_input == -1 || bp_input < bp_ref) && snp_input < input.GetL()){
      input.ReadSNP(sequence_input, bp_input);
      snp_input++;
      if(snp_input == input.GetL()) break;
    }
    while(bp_ref > mut.info[snp].pos){
      snp++;
      if(snp == mut.info.size()) break;
    }
    if(snp == mut.info.size()) break;

    if(bp_ref == mut.info[snp].pos){

      assert(bp_ref == bp_input);
      int bin_index = std::max(0, (int)std::round(log(10*age[snp])*C));
      assert(bin_index < num_age_bins);
      bool used = false;
      //calculate contribution to MLE if shared and non-shared
      for(int i = 0; i < input.GetN(); i++){
        for(int j = 0; j < ref.GetN(); j++){
          //int i = 0, j = 15;
          if(mut.info[snp].age_end > 0 && age[snp] > 0 && age[snp] > 0 && mut.info[snp].age_end < tmrca[snp]){
            if(sequence_ref[j] == '1'){

              if(!used){
                used = true;
                num_used_snps++;
              }
              if(bp_ref == bp_input && sequence_input[i] == '1'){
                age_shared_count[j][bin_index]++;
              }else{
                age_notshared_count[j][bin_index]++;
              }

            }
          }
        }
      }

    }else{
      //assert(1);
    }
  }
  ref.CloseFile();
  input.CloseFile();

  //std::cerr << num_used_snps << " " << bp_ref << std::endl;

  //////////////////

  int max_iter = 1000;
  int perc = -1;
  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      std::cerr << "[" << perc << "%]\r"; 
    }

    //double start_time = time(NULL);
    //clock_t begin = clock();
    //coal_EM_simplified EM(epochs, coal_rates[0]);
    std::vector<coal_EM_simplified> EM;
    EM.reserve(data.N);
    for(int i = 0; i < data.N; i++){
      EM.push_back(coal_EM_simplified(epochs, coal_rates[i]));
    }

    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;

    //start_time = time(NULL);
    //begin = clock(); 
    int count = 0;  
    std::vector<double> num(num_epochs,0.0), denom(num_epochs,0.0);

    for(int bin = 0; bin < num_age_bins; bin++){
      for(int j = 0; j < data.N; j++){
        if(age_shared_count[j][bin] > 0){ 
          count = age_shared_count[j][bin];
          std::fill(num.begin(), num.end(), 0.0);
          std::fill(denom.begin(), denom.end(), 0.0);
          double logl = EM[j].EM_shared(age_bin[bin], num, denom);
          //assert(logl <= 0.0);
          log_likelihood += count * logl;
          for(int e = 0; e < num_epochs-1; e++){
            coal_rates_num[j][e]   += count * num[e];
            coal_rates_denom[j][e] += count * denom[e];
          }
        }
        if(age_notshared_count[j][bin] > 0){ 
          count = age_notshared_count[j][bin];
          std::fill(num.begin(), num.end(), 0.0);
          std::fill(denom.begin(), denom.end(), 0.0);
          double logl = EM[j].EM_notshared(age_bin[bin], num, denom);
          //assert(logl <= 0.0);
          log_likelihood += count * logl;
          for(int e = 0; e < num_epochs-1; e++){
            coal_rates_num[j][e]   += count * num[e];
            coal_rates_denom[j][e] += count * denom[e];
          }
        }
      }
    }
    //end = clock();
    //end_time = time(NULL);
    //double elapsed_secs2 = double(end - begin) / CLOCKS_PER_SEC;
    //std::cerr << elapsed_secs1 << " " << elapsed_secs2 << std::endl;

    for(int i = 0; i < data.N; i++){
      for(int e = 0; e < num_epochs; e++){
        if(coal_rates_num[i][e] == 0){
          if(e > 0){
            coal_rates[i][e] = coal_rates[i][e-1];
          }else{
            assert(1);
            coal_rates[i][e] = 0;
          }
        }else if(coal_rates_denom[i][e] == 0){
          assert(1);
        }else{
          coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
        }
      }
    }

    //std::cerr << std::endl;
    os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << "\n";

    for(int i = 0; i < data.N; i++){
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }

  }
  os_log.close();

  //output
  std::ofstream os(options["output"].as<std::string>() + ".coal");
  for(int i = 0; i < data.N; i++){
    os << i << " ";
  }
  os << std::endl;

  for(int e = 0; e < num_epochs; e++){
    os << epochs[e] << " ";
  }
  os << std::endl;

  for(int i = 0; i < data.N; i++){
    os << "0 " << i << " ";
    for(int e = 0; e < num_epochs; e++){
      os << coal_rates[i][e] << " ";
    }
    os << std::endl;
  }

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
mut_fast_simplified_all(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, haps, sample, input, output. Optional: num_bins, coal" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Calculating coalescence rates for (ancient) sample.." << std::endl;

  ////////////////////////////////////////

  int regularise = 0;
  if(options.count("regularise") > 0){ 
    regularise = options["regularise"].as<int>();
  }

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
      assert(epochs[e] > epochs[e-1]);
    }

  }else if(1){
    num_epochs = 20;
    if(options.count("num_bins") > 0){
      num_epochs = options["num_bins"].as<int>();
    }
    float years_per_gen = 28.0;
    if(options.count("years_per_gen")){
      years_per_gen = options["years_per_gen"].as<float>();
    }
    num_epochs++; 
    epochs.resize(num_epochs);

    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }else if(0){
    num_epochs = 6;
    epochs.resize(num_epochs);

    epochs[0] = 0;
    epochs[1] = 1e3/28;
    epochs[2] = 1e4/28;
    epochs[3] = 1e5/28;
    epochs[4] = 0.99e7/28;
    epochs[5] = 1e8;

  }else{
    num_epochs = 4;
    epochs.resize(num_epochs);
    epochs[0] = 0;
    epochs[1] = 1e4/28;
    epochs[2] = 1e7/28;
    epochs[3] = 1e8/28;
  }

  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)
  //For now, assume both are haps/sample and span same positions

  //haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  //Data data(ref_N.GetN(), 1);
  //ref_N.CloseFile();
  std::vector<double> coal_rates(num_epochs), coal_rates_num(num_epochs), coal_rates_denom(num_epochs);

  if(options.count("coal") > 0){
    double dummy;
    is >> dummy >> dummy;
    for(int e = 0; e < num_epochs; e++){
      is >> coal_rates[e];
      std::cerr << coal_rates[e] << " "; 
    }
    std::cerr << std::endl;
  }else{
    std::fill(coal_rates.begin(), coal_rates.end(), 1.0/100000.0);
    std::fill(coal_rates_num.begin(), coal_rates_num.end(), 0.0);
    std::fill(coal_rates_denom.begin(), coal_rates_denom.end(), 0.0);
  }


  ////////////////////////////////////////  

  std::ofstream os_log(options["output"].as<std::string>() + ".log");

  //for every input & ref haplotype pair, want to know where they share and not share
  //for ref haplotype, record age of SNPs that share & age of SNPs that don't share
  //const*N_ref*N_input
  //could even bin this by discretised age if wanted

  //formula for converting age to int: log(age)*C, Ne = haploid population size
  double C = 1e2;
  int num_age_bins = ((int) (log(1e8) * C)+1);
  //std::cerr << num_age_bins << std::endl;
  std::vector<double> age_bin(num_age_bins, 0.0);
  int bin = 0;
  for(std::vector<double>::iterator it_age_bin = std::next(age_bin.begin(),1); it_age_bin != age_bin.end(); it_age_bin++){
    *it_age_bin = exp(bin/C)/10.0;
    bin++;
  }

  std::vector<double> age_shared_count(num_age_bins), age_notshared_count(num_age_bins);
  std::fill(age_shared_count.begin(), age_shared_count.end(), 0.0);
  std::fill(age_notshared_count.begin(), age_notshared_count.end(), 0.0);

  ///////////

  std::cerr << "Parsing input files" << std::endl;

  /////////////////////////

  int num_used_snps = 0;
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
    std::cerr << chromosomes.size() << std::endl;

  }else{
    chromosomes.resize(0);
    exit(1);
  }

  for(int chr = 0; chr < chromosomes.size(); chr++){

    std::cerr << "CHR: " << chr << std::endl;

    std::string filename_anc      = options["anc"].as<std::string>() + "_chr" + chromosomes[chr] + ".anc";
    std::string filename_mut      = options["mut"].as<std::string>() + "_chr" + chromosomes[chr] + ".mut";
    std::string filename_haps     = options["haps"].as<std::string>() + "_chr" + chromosomes[chr] + ".haps.gz";
    std::string filename_sample   = options["sample"].as<std::string>() + "_chr" + chromosomes[chr] + ".sample.gz";
    std::string file_input_haps   = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".haps.gz";
    std::string file_input_sample = options["input"].as<std::string>() + "_chr" + chromosomes[chr] + ".sample.gz";

    Mutations mut;
    mut.Read(filename_mut);
    haps ref(filename_haps.c_str(), filename_sample.c_str());
    haps input(file_input_haps.c_str(), file_input_sample.c_str());

    //assign an actual date to each SNP
    std::uniform_real_distribution<double> dist_unif(0,1);   
    std::mt19937 rng;
    std::vector<double> age(mut.info.size(), 0.0);
    std::vector<double>::iterator it_age = age.begin();
    int seed = 10;
    rng.seed(seed);
    for(Muts::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++ ){
      *it_age = dist_unif(rng) * ((*it_mut).age_end - (*it_mut).age_begin) + (*it_mut).age_begin;
      //*it_age = ((*it_mut).age_end + (*it_mut).age_begin)/2.0;
      it_age++;
    }

    int bp_ref, bp_input = -1;
    int snp = 0, snp_input = 0, snp_ref = 0;
    std::vector<char> sequence_ref(ref.GetN()), sequence_input(input.GetN());
    //iterating over snps in reference
    for(; snp_ref < ref.GetL();){

      for(int k = 0; k < 1; k++){
        ref.ReadSNP(sequence_ref, bp_ref);
        snp_ref++;
        if(snp_ref == ref.GetL()) break;
      }

      int DAF = 0;
      for(std::vector<char>::iterator it_seq = sequence_ref.begin(); it_seq != sequence_ref.end(); it_seq++){
        DAF += (*it_seq == '1');
      }

      //std::cerr << bp_ref << std::endl;
      snp_ref++;
      while((bp_input == -1 || bp_input < bp_ref) && snp_input < input.GetL()){
        input.ReadSNP(sequence_input, bp_input);
        snp_input++;
        if(snp_input == input.GetL()) break;
      }
      while(bp_ref > mut.info[snp].pos){
        snp++;
        if(snp == mut.info.size()) break;
      }
      if(snp == mut.info.size()) break;

      if(bp_ref == mut.info[snp].pos){

        assert(bp_ref == bp_input);
        int bin_index = std::max(0, (int)std::round(log(10*age[snp])*C)+1);
        assert(bin_index < num_age_bins);
        bool used = false;
        //calculate contribution to MLE if shared and non-shared
        for(int i = 0; i < input.GetN(); i++){
          for(int j = 0; j < ref.GetN(); j++){
            if(mut.info[snp].age_end > 0 && mut.info[snp].age_end > mut.info[snp].age_begin && mut.info[snp].flipped == 0 && DAF > 0){

              if(sequence_ref[j] == '1'){

                if(!used){
                  used = true;
                  num_used_snps++;
                }

                if(1){
                  if(bp_ref == bp_input && sequence_input[i] == '1'){
                    age_shared_count[bin_index] += 1;
                  }else{
                    age_notshared_count[bin_index] += 1;
                  }                
                }else{
                  int max_k = 100;
                  for(int k = 0; k < max_k; k++){

                    //double snp_age = dist_unif(rng) * (mut.info[snp].age_end - mut.info[snp].age_begin) + mut.info[snp].age_begin;
                    double snp_age = (k+1.0)/((double)max_k) * (mut.info[snp].age_end - mut.info[snp].age_begin) + mut.info[snp].age_begin;
                    bin_index = std::max(0, (int)std::round(log(10*snp_age)*C)+1);
                    double weight = 1.0/max_k;
                    //weight =  mut.info[snp].age_end - mut.info[snp].age_begin;
                    if(bp_ref == bp_input && sequence_input[i] == '1'){
                      age_shared_count[bin_index] += weight;
                    }else{
                      age_notshared_count[bin_index] += weight;
                    }

                  }
                }

              }

            }
          }
        }

      }
    }
    ref.CloseFile();
    input.CloseFile();
  }

  std::cerr << num_used_snps << std::endl;

  double gamma = 1;
  std::vector<double> gradient(num_epochs, 0.0), gradient_prev(num_epochs, 0.0), coal_rates_prev(num_epochs, 0.0), coal_rates_nesterov(num_epochs, 0.0);
  coal_rates_nesterov = coal_rates;

  int max_iter = 10000;
  int perc = -1;
  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    /*
       std::cerr << iter << " ";
       for(int e = 0; e < num_epochs; e++){
       std::cerr << coal_rates[e] << " ";
       }
       std::cerr << std::endl;
       */
    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      std::cerr << "[" << perc << "%]\r"; 
    }

    //double start_time = time(NULL);
    //clock_t begin = clock();
    //coal_EM_simplified EM(epochs, coal_rates[0]);
    coal_EM_simplified EM(epochs, coal_rates_nesterov);

    //clock_t end = clock();
    //double end_time = time(NULL);
    //double elapsed_secs1 = double(end - begin) / CLOCKS_PER_SEC;

    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;

    //start_time = time(NULL);
    //begin = clock(); 
    double count = 0;  
    std::vector<double> num(num_epochs,0.0), denom(num_epochs,0.0);
    std::vector<double> num_diff(num_epochs,0.0), denom_diff(num_epochs,0.0);
    double logl_diff = 0.0;

    for(int bin = 1; bin < num_age_bins; bin++){
      if(age_shared_count[bin] > 0){ 
        count = age_shared_count[bin];

        std::fill(num.begin(), num.end(), 0.0);
        std::fill(denom.begin(), denom.end(), 0.0);
        double logl;
        logl = EM.EM_shared(age_bin[bin], num, denom);
        assert(logl <= 0.0);
        log_likelihood += count * logl;
        for(int e = 0; e < num_epochs-1; e++){
          coal_rates_num[e]   += count * num[e];
          coal_rates_denom[e] += count * denom[e];
        }
      }
      if(age_notshared_count[bin] > 0){ 
        count = age_notshared_count[bin];
        std::fill(num.begin(), num.end(), 0.0);
        std::fill(denom.begin(), denom.end(), 0.0);
        double logl;
        logl = EM.EM_notshared(age_bin[bin], num, denom);
        assert(logl <= 0.0);
        log_likelihood += count * logl;
        for(int e = 0; e < num_epochs-1; e++){
          //num_diff[e] += count * (num[e]-num2[e]);
          //denom_diff[e] += count * (denom[e]-denom2[e]);
          coal_rates_num[e]   += count * num[e];
          coal_rates_denom[e] += count * denom[e];
        }
      }
    }
    //end = clock();
    //end_time = time(NULL);
    //double elapsed_secs2 = double(end - begin) / CLOCKS_PER_SEC;
    //std::cerr << elapsed_secs1 << " " << elapsed_secs2 << std::endl;

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
      gamma = 1e-11/num_used_snps;
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
        if(coal_rates[e] < 1e-7){
          coal_rates[e] = 1e-7;
        }
        if(coal_rates_nesterov[e] < 1e-7){
          coal_rates_nesterov[e] = 1e-7;
        }
        assert(coal_rates[e] >= 0.0);
        assert(coal_rates_nesterov[e] >= 0.0);

      }

      if(!is_EM){
        //gradient_prev[e]   = gradient[e];
        gradient_prev[e] = 0.9*gradient_prev[e] + gamma * gradient[e];
      }

    }

    //std::cerr << std::endl;
    os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << "\n";

    std::fill(coal_rates_num.begin(), coal_rates_num.end(), 0.0);
    std::fill(coal_rates_denom.begin(), coal_rates_denom.end(), 0.0);

  }
  os_log.close();

  //output
  std::ofstream os(options["output"].as<std::string>() + ".coal");

  os << "0\n";
  for(int e = 0; e < num_epochs; e++){
    os << epochs[e] << " ";
  }
  os << std::endl;

  os << "0 0" << " ";
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
  if(!options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mut, haps, sample, input, output." << std::endl;
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
  int root = 2*N-2;

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
  float num_bases_snp_persists = 0.0;
  std::vector<float> coords;
  double tmrca = 0.0;

  Mutations mut_out;
  mut_out.Read(options["input"].as<std::string>() + ".mut");

  Mutations mut_combined;
  mut_combined.info.resize(L);

  int L_ref = ancmut.NumSnps();
  int L_out = mut_out.info.size();

  int snp_ref = 0, snp_out = 0, snp_hap = 0;
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
    if(snp_out < L_out){
      while(mut_out.info[snp_out].pos < bp){
        snp_out++;
        if(snp_out == L_out) break;
      }
    }

    if((*it_mut).pos == bp && DAF > 0 && DAF < N){
      //if segregating, copy over
      if(tree_count < (*it_mut).tree){
        tree_count = (*it_mut).tree;
        mtr.tree.GetCoordinates(coords);
        tmrca = coords[root];
      }
      mut_combined.info[snp_hap] = (*it_mut);
      if(mut_combined.info[snp_hap].branch.size() == 1){
        assert( *mut_combined.info[snp_hap].branch.begin() <= root );
      }
    }else{
      //otherwise copy from trees with outgroup
      if(0){
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
      }else{      
        if(DAF == N){
          if(tmrca <= outgroup_age){
            mut_combined.info[snp_hap].branch.resize(1);
            mut_combined.info[snp_hap].branch[0] = root;
            mut_combined.info[snp_hap].age_begin = tmrca;
            mut_combined.info[snp_hap].age_end   = outgroup_age;
            assert(mut_combined.info[snp_hap].age_begin <= mut_combined.info[snp_hap].age_end);
          }
        }else{
          mut_combined.info[snp_hap].age_begin = 0;
          mut_combined.info[snp_hap].age_end = 0;
          mut_combined.info[snp_hap].branch.clear();
        }
      }
    }

    mut_combined.info[snp_hap].tree = tree_count;
    mut_combined.info[snp_hap].pos = bp;
    if(snp_hap > 0) mut_combined.info[snp_hap-1].dist = bp - mut_combined.info[snp_hap-1].pos; 
    mut_combined.info[snp_hap].snp_id = snp_hap;	
  }

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


