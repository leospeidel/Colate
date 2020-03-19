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

#include "aDNA_EM.hpp"

void
aDNA(cxxopts::Options& options){

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

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());

  /////////////////////////////////
  //get TMRCA at each SNP
  std::vector<float> tmrca(mut.info.size(), 1e8/28.0);
  std::vector<float>::iterator it_tmrca = tmrca.begin();
  std::vector<double> age_begin(mut.info.size(), 0.0), age_end(mut.info.size(), 0.0);

  if(0){
    MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
    Muts::iterator it_mut; //iterator for mut file
    AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
    std::vector<double>::iterator it_age_begin = age_begin.begin(), it_age_end = age_end.begin();

    float num_bases_tree_persists = 0.0;
    int root = 2 * ancmut.NumTips() - 2;
    std::vector<float> coords(2*ancmut.NumTips() - 1);
    while(num_bases_tree_persists >= 0.0){
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

      if(num_bases_tree_persists >= 0.0){
        mtr.tree.GetCoordinates(coords);
        int tree_index = (*it_mut).tree;
        while((*it_mut).tree == tree_index){
          *it_tmrca = coords[root];
          if((*it_mut).branch.size() == 1 && (*it_mut).flipped == 0){
            *it_age_begin = coords[*(*it_mut).branch.begin()];
            *it_age_end   = coords[(*mtr.tree.nodes[*(*it_mut).branch.begin()].parent).label];
          }

          it_age_begin++;
          it_age_end++;

          it_tmrca++;
          it_mut++;
        }
      }

    }
    std::ofstream os_tmrca(options["input"].as<std::string>() + ".tmrca");
    for(it_tmrca = tmrca.begin(); it_tmrca != tmrca.end(); it_tmrca++){
      os_tmrca << *it_tmrca << "\n";
    }
    os_tmrca.close();
  }

  if(0){
    std::ifstream is_tmrca(options["input"].as<std::string>() + ".tmrca");
    for(it_tmrca = tmrca.begin(); it_tmrca != tmrca.end(); it_tmrca++){
      is_tmrca >> *it_tmrca;
    }
    is_tmrca.close();
  }

  std::vector<float> tmrca2(mut.info.size(), 1e8/28.0);
  std::vector<float>::iterator it_tmrca2 = tmrca2.begin();
  if(0){
    std::ifstream is_tmrca2(options["input"].as<std::string>() + ".tmrca");
    for(it_tmrca2 = tmrca2.begin(); it_tmrca2 != tmrca2.end(); it_tmrca2++){
      is_tmrca2 >> *it_tmrca2;
    }
    is_tmrca2.close();
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
    num_epochs = 10;
    if(options.count("num_bins") > 0){
      num_epochs = options["num_bins"].as<int>();
    }
    float years_per_gen = 28.0;
    if(options.count("years_per_gen")){
      years_per_gen = options["years_per_gen"].as<float>();
    }
    num_epochs++; 
    epochs.resize(num_epochs);

    if(1){
      epochs[0] = 0.0;
      epochs[1] = 1e3/years_per_gen;
      float log_10 = std::log(10);
      for(int e = 2; e < num_epochs-1; e++){
        epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
      }
      epochs[num_epochs-1] = 1e8/years_per_gen;
    }else{
      num_epochs = 7;
      epochs.resize(num_epochs);
      epochs[0] = 0.0;
      epochs[1] = 5e4/28.0;
      epochs[2] = 1e5/28.0;
      epochs[3] = 5e5/28.0;
      epochs[4] = 1e6/28.0;
      epochs[5] = 5e6/28.0;
      epochs[6] = 1e8/28.0;
    }
  }


  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)
  //For now, assume both are haps/sample and span same positions
  //TODO: make this more flexible

  haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  Data data(ref_N.GetN(), 1);
  ref_N.CloseFile();
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
    double initial_coal_rate = 1.0/40000.0;
    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
      /*
         coal_rates[i][0] = 1.43472e-05;
         coal_rates[i][1] = 1.34137e-05;
         coal_rates[i][2] = 1.65727e-05;
         coal_rates[i][3] = 2.05526e-05;
         coal_rates[i][4] = 5.98515e-05;
         coal_rates[i][5] = 0.000118802;
         */
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }
  }

  ////////////////////////////////////////  

  int N = data.N;
  int L = mut.info.size();
  int max_iter = 1000;
  int perc = -1;
  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      std::cerr << "[" << perc << "%]\r";
    }

    std::vector<aDNA_EM> EM;
    EM.reserve(N);
    for(int i = 0; i < N; i++){
      EM.push_back(aDNA_EM(epochs, coal_rates[i]));
    }

    haps ref(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
    haps input((options["input"].as<std::string>() + ".haps.gz").c_str(), (options["input"].as<std::string>() + ".sample.gz").c_str());

    int bp_ref, bp_input = -1;
    std::vector<char> sequence_ref(data.N), sequence_input(input.GetN());
    std::vector<float> num(num_epochs, 0.0), denom(num_epochs, 0.0); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1

    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;
    int snp = 0, snp_input = 0, snp_ref = 0;

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
        if(snp == L) break;
      }
      if(snp == L) break;

      if(bp_ref == mut.info[snp].pos){

        //calculate contribution to MLE if shared and non-shared
        for(int i = 0; i < N; i++){
          if(mut.info[snp].age_end > 0){
            if(sequence_ref[i] == '1'){

              //int ep = 0;
              //while(epochs[ep] < mut.info[snp].age_end){
              //  ep++;
              //  if(ep == epochs.size()) break;
              //}

              //if(ep+1 < epochs.size()){
              //if(tmrca2[snp] > epochs[ep+1]){

              if(bp_ref == bp_input && sequence_input[0] == '1'){
                log_likelihood += EM[i].EM_shared(mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], num, denom);
                //EM[i].EM_shared(age_begin[snp], age_end[snp], tmrca[snp], num, denom);
                assert(!std::isnan(denom[0]));
              }else{
                log_likelihood += EM[i].EM_notshared(mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], num, denom);
                //EM[i].EM_notshared(age_begin[snp], age_end[snp], tmrca[snp], num, denom);
                assert(!std::isnan(denom[0]));
              }


              for(int e = 0; e < num_epochs; e++){
                coal_rates_num[i][e]   += num[e];
                coal_rates_denom[i][e] += denom[e];
              }

              /*
                 for(int e = 0; e < num_epochs; e++){
                 coal_rates_num[0][e]   += num[e]/N;
                 coal_rates_denom[0][e] += denom[e]/N;
                 }
                 */

            }
            //}

            //}
          }
        }

      }

    }

    for(int i = 0; i < N; i++){
      for(int e = 0; e < num_epochs; e++){
        if(coal_rates_num[i][e] == 0){
          coal_rates[i][e] = 0;
        }else if(coal_rates_denom[i][e] == 0){
        }else{
          coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
        }
      }
    }

    /*   
         for(int i = 0; i < N; i++){
         for(int e = 0; e < num_epochs; e++){
         if(coal_rates_num[0][e] == 0){
         coal_rates[0][e] = 0;
         }else if(coal_rates_denom[0][e] == 0){
         assert(1);
         }else{
         coal_rates[i][e] = coal_rates_num[0][e]/coal_rates_denom[0][e];
         }
         }
         }
         */

    ref.CloseFile();
    input.CloseFile();

    std::cerr << log_likelihood << " " << log_likelihood - prev_log_likelihood << std::endl;

    /*    
          for(int e = 0; e < num_epochs; e++){
          std::cerr << coal_rates[0][e] << " ";
          }
          std::cerr << std::endl;
          */

    if(log_likelihood - prev_log_likelihood < 0.0) break;
    /*
       for(int i = 0; i < N; i++){
       for(int e = 0; e < num_epochs; e++){
       std::cerr << coal_rates[i][e] << " ";
       }
       std::cerr << std::endl;
       }
       */

    for(int i = 0; i < N; i++){
      std::fill(coal_rates_num[i].begin(), coal_rates_num[i].end(), 0.0);
      std::fill(coal_rates_denom[i].begin(), coal_rates_denom[i].end(), 0.0);
    }

    std::ofstream os(options["output"].as<std::string>());

    for(int i = 0; i < N; i++){
      os << i << " ";
    }
    os << std::endl;

    for(int e = 0; e < num_epochs; e++){
      os << epochs[e] << " ";
    }
    os << std::endl;

    for(int i = 0; i < N; i++){
      os << "0 " << i << " ";
      for(int e = 0; e < num_epochs; e++){
        os << coal_rates[i][e] << " ";
      }
      os << std::endl;
    }

    os.close();

  }

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
aDNA_simplified(cxxopts::Options& options){

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
  int seed = 1;
  rng.seed(seed);
  for(Muts::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++ ){
    *it_age = dist_unif(rng) * ((*it_mut).age_end - (*it_mut).age_begin) + (*it_mut).age_begin;
    it_age++;
  }

  /////////////////////////////////
  //get TMRCA at each SNP
  std::vector<float> tmrca(mut.info.size(), 1e8/28.0);
  std::vector<float>::iterator it_tmrca = tmrca.begin();
  std::vector<double> age_begin(mut.info.size(), 0.0), age_end(mut.info.size(), 0.0);

  bool correct = 0;
  if(options.count("correction") > 0){ 
    correct = options["correction"].as<int>();
  }

  if(correct){
    MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
    Muts::iterator it_mut; //iterator for mut file
    AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
    std::vector<double>::iterator it_age_begin = age_begin.begin(), it_age_end = age_end.begin();

    float num_bases_tree_persists = 0.0;
    int root = 2 * ancmut.NumTips() - 2;
    std::vector<float> coords(2*ancmut.NumTips() - 1);
    while(num_bases_tree_persists >= 0.0){
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

      if(num_bases_tree_persists >= 0.0){
        mtr.tree.GetCoordinates(coords);
        int tree_index = (*it_mut).tree;
        while((*it_mut).tree == tree_index){
          *it_tmrca = coords[root];
          if((*it_mut).branch.size() == 1 && (*it_mut).flipped == 0){
            *it_age_begin = coords[*(*it_mut).branch.begin()];
            *it_age_end   = coords[(*mtr.tree.nodes[*(*it_mut).branch.begin()].parent).label];
          }

          it_age_begin++;
          it_age_end++;

          it_tmrca++;
          it_mut++;
        }
      }

    }
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

  //haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  //Data data(ref_N.GetN(), 1);
  //ref_N.CloseFile();
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
    double initial_coal_rate = 1.0/40000.0;
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

  haps ref(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  haps input((options["input"].as<std::string>() + ".haps.gz").c_str(), (options["input"].as<std::string>() + ".sample.gz").c_str());

  std::vector<std::vector<std::vector<double>>> age_shared(input.GetN());
  std::vector<std::vector<std::vector<double>>> age_notshared(input.GetN());
  for(int i = 0; i < input.GetN(); i++){
    age_shared[i].resize(ref.GetN());
    age_notshared[i].resize(ref.GetN());
    for(int j = 0; j < ref.GetN(); j++){
      age_shared[i][j].reserve(ref.GetL());
      age_notshared[i][j].reserve(ref.GetL());
    }
  }

  int bp_ref, bp_input = -1;
  int snp = 0, snp_input = 0, snp_ref = 0;
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
      if(snp == ref.GetL()) break;
    }
    if(snp == ref.GetL()) break;

    if(bp_ref == mut.info[snp].pos){

      //calculate contribution to MLE if shared and non-shared
      for(int i = 0; i < input.GetN(); i++){
        for(int j = 0; j < ref.GetN(); j++){
          if(mut.info[snp].age_end > 0 && mut.info[snp].age_end < tmrca[snp]){
            if(sequence_ref[j] == '1'){

              if(bp_ref == bp_input && sequence_input[i] == '1'){
                age_shared[i][j].push_back(age[snp]);
              }else{
                age_notshared[i][j].push_back(age[snp]);
              }

            }
          }
        }
      }

    }
  }

  //////////////////

  int max_iter = 1000;
  int perc = -1;
  double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
  for(int iter = 0; iter < max_iter; iter++){

    if( (int) (((double)iter)/max_iter * 100.0) > perc ){
      perc = (int) (((double)iter)/max_iter * 100.0);
      std::cerr << "[" << perc << "%]\r";

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
    }

    double start_time = time(NULL);
    clock_t begin = clock();
    std::vector<aDNA_EM_simplified> EM;
    EM.reserve(data.N);
    for(int i = 0; i < data.N; i++){
      EM.push_back(aDNA_EM_simplified(epochs, coal_rates[i]));
    }
    clock_t end = clock();
    double end_time = time(NULL);
    double elapsed_secs1 = double(end - begin) / CLOCKS_PER_SEC;

    prev_log_likelihood = log_likelihood;
    log_likelihood = 0.0;

    start_time = time(NULL);
    begin = clock();
    for(int i = 0; i < input.GetN(); i++){
      for(int j = 0; j < ref.GetN(); j++){
        for(std::vector<double>::iterator it_age = age_shared[i][j].begin(); it_age != age_shared[i][j].end(); it_age++){
          log_likelihood += EM[0].EM_shared(*it_age, coal_rates_num[0], coal_rates_denom[0]);
        }
        for(std::vector<double>::iterator it_age = age_notshared[i][j].begin(); it_age != age_notshared[i][j].end(); it_age++){
          log_likelihood += EM[0].EM_notshared(*it_age, coal_rates_num[0], coal_rates_denom[0]); 
        }
      }
    }     
    end = clock();
    end_time = time(NULL);
    double elapsed_secs2 = double(end - begin) / CLOCKS_PER_SEC;
    //std::cerr << elapsed_secs1 << " " << elapsed_secs2 << std::endl;
    /*
       for(int i = 0; i < N; i++){
       for(int e = 0; e < num_epochs; e++){
       if(coal_rates_num[i][e] == 0){
       coal_rates[i][e] = 0;
       }else if(coal_rates_denom[i][e] == 0){
       }else{
       coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
       }
       }
       }
       */

    for(int i = 0; i < data.N; i++){
      for(int e = 0; e < num_epochs; e++){
        if(coal_rates_num[i][e] == 0){
          coal_rates[i][e] = 0;
        }else if(coal_rates_denom[i][e] == 0){
          assert(1);
        }else{
          coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
        }
      }
    }

    os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << std::endl;

    if(log_likelihood - prev_log_likelihood < 0.0) break;

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
aDNA_fast_simplified(cxxopts::Options& options){

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
          if(it_mut == mut.info.end()) break;
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
	//for(;snp < 217443;){
  //for(;snp < 300000;){
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
      //for(int i = 0; i < input.GetN(); i++){
			int i = 0;
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
      //}

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
    //aDNA_EM_simplified EM(epochs, coal_rates[0]);
    std::vector<aDNA_EM_simplified> EM;
    EM.reserve(data.N);
    for(int i = 0; i < data.N; i++){
      EM.push_back(aDNA_EM_simplified(epochs, coal_rates[i]));
    }

    //clock_t end = clock();
    //double end_time = time(NULL);
    //double elapsed_secs1 = double(end - begin) / CLOCKS_PER_SEC;

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

    //if(log_likelihood - prev_log_likelihood < 0.0) break;

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
aDNA_fast_simplified_all(cxxopts::Options& options){

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

  /*
	if(0){
		//recalculate ages
		MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
		Muts::iterator it_mut; //iterator for mut file
		AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());

		it_age = age.begin();
		float num_bases_tree_persists = 0.0;
		int root = 2 * ancmut.NumTips() - 2;
		std::vector<float> coords(2*ancmut.NumTips() - 1);
		int tree_count = 0;
		int snp = 0;
		while(num_bases_tree_persists >= 0.0){
			num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
			if(num_bases_tree_persists >= 0.0 && (*it_mut).tree == tree_count){
				mtr.tree.GetCoordinates(coords);
				int tree_index = (*it_mut).tree;
				while((*it_mut).tree == tree_index){
					if((*it_mut).branch.size() == 1){
					  int node = *(*it_mut).branch.begin();
					  *it_age = dist_unif(rng) * ( coords[(*mtr.tree.nodes[node].parent).label] - coords[node] ) + coords[node];
					}else{
            *it_age = 0.0;
					}
					it_age++;
					it_mut++;
					snp++;
					if(it_mut == mut.info.end()) break;
				}
			}
			tree_count++;
		}
		ancmut.CloseFiles();
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
					if(it_mut == mut.info.end()) break;
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
  */

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

	//haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
	//Data data(ref_N.GetN(), 1);
	//ref_N.CloseFile();
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
		double initial_coal_rate = 1.0/40000.0;
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

	/////////////////////////

	Mutations mut;

  for(int chr = 1; chr <= 30; chr++){
	mut.Read("./chr_" + std::to_string(chr) + "/" + options["mut"].as<std::string>());
	haps ref(("./chr_" + std::to_string(chr) + "/" +options["haps"].as<std::string>()).c_str(), ("./chr_" + std::to_string(chr) + "/" + options["sample"].as<std::string>()).c_str());
	haps input(("./chr_" + std::to_string(chr) + "/" + options["input"].as<std::string>() + ".haps.gz").c_str(), ("./chr_" + std::to_string(chr) + "/" + options["input"].as<std::string>() + ".sample.gz").c_str());


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

	int bp_ref, bp_input = -1;
	int snp = 0, snp_input = 0, snp_ref = 0;
	std::vector<char> sequence_ref(ref.GetN()), sequence_input(input.GetN());
	//iterating over snps in reference
	int num_used_snps = 0;
	for(; snp_ref < ref.GetL();){
	
		for(int k = 0; k < 1; k++){
		  ref.ReadSNP(sequence_ref, bp_ref);
		  snp_ref++;
			if(snp_ref == ref.GetL()) break;
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
			int bin_index = std::max(0, (int)std::round(log(10*age[snp])*C));
			assert(bin_index < num_age_bins);
			bool used = false;
			//calculate contribution to MLE if shared and non-shared
			//for(int i = 0; i < input.GetN(); i++){
			int i = 0;
			for(int j = 0; j < ref.GetN(); j++){
			//int j = 0;
				//int i = 0, j = 15;
				//if(mut.info[snp].age_end > 0 && age[snp] > epochs[0] && age[snp] < epochs[1] && mut.info[snp].age_end < tmrca[snp]){
				if(mut.info[snp].age_end > 0){
					if(sequence_ref[j] == '1'){

						if(!used){
							used = true;
							num_used_snps++;
						}
						if(bp_ref == bp_input && sequence_input[i] == '1'){
							age_shared_count[0][bin_index]++;
						}else{
							age_notshared_count[0][bin_index]++;
						}

					}
				}
			}
			//}

		}else{
			//assert(1);
		}
	}
	ref.CloseFile();
	input.CloseFile();
  }

	//std::cerr << num_used_snps << " " << bp_ref << std::endl;

	//////////////////

	/*
	aDNA_EM_simplified EM_tmp(epochs, coal_rates[0]);
	for(int i = 0; i < num_epochs-2; i++){
		std::vector<double> num(num_epochs,0.0), denom(num_epochs,0.0);
	  int bin_min = std::max(0, (int)std::round(log(10*epochs[i])*C)), bin_max = std::max(0, (int)std::round(log(10*epochs[i+1])*C));
		int count1 = 0, count2 = 0;
		std::cerr << "epoch: " << epochs[i] << std::endl;
		for(int bin = bin_min; bin < bin_max; bin++){
			std::fill(num.begin(), num.end(), 0.0);
			std::fill(denom.begin(), denom.end(), 0.0);
			double logl = EM_tmp.EM_notshared(age_bin[bin], num, denom);
			std::cerr << bin << " " << num[0] << " " << denom[0] << " " << epochs[1] - epochs[0] << " " << num[0]*age_shared_count[0][bin] << " " << denom[0]*age_notshared_count[0][bin] << std::endl;
      count1 += age_shared_count[0][bin];
			count2 += age_notshared_count[0][bin];
		}
	  //std::cerr << count1 << " " << count2 << std::endl;
  }
	std::cerr << std::endl;
	*/

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
		//aDNA_EM_simplified EM(epochs, coal_rates[0]);
		std::vector<aDNA_EM_simplified> EM;
		EM.reserve(data.N);
		for(int i = 0; i < data.N; i++){
			EM.push_back(aDNA_EM_simplified(epochs, coal_rates[i]));
		}

		//clock_t end = clock();
		//double end_time = time(NULL);
		//double elapsed_secs1 = double(end - begin) / CLOCKS_PER_SEC;

		prev_log_likelihood = log_likelihood;
		log_likelihood = 0.0;

		//start_time = time(NULL);
		//begin = clock(); 
		int count = 0;  
		std::vector<double> num(num_epochs,0.0), denom(num_epochs,0.0);
		std::vector<double> num_diff(num_epochs,0.0), denom_diff(num_epochs,0.0);
		double logl_diff = 0.0;

		for(int bin = 0; bin < num_age_bins; bin++){
			for(int j = 0; j < data.N; j++){
				if(age_shared_count[j][bin] > 0){ 
					count = age_shared_count[j][bin];
					
					std::fill(num.begin(), num.end(), 0.0);
					std::fill(denom.begin(), denom.end(), 0.0);
					double logl;
					std::vector<double> num2(num_epochs,0.0), denom2(num_epochs,0.0);
					logl = EM[j].EM_shared_exact(age_bin[bin], num2, denom2);
					logl_diff += logl - EM[j].EM_shared(age_bin[bin], num, denom);
					assert(logl <= 0.0);
					log_likelihood += count * logl;
					for(int e = 0; e < num_epochs-1; e++){
						num_diff[e]   += count * (num[e]-num2[e]);
						denom_diff[e] += count * (denom[e]-denom2[e]);
						coal_rates_num[j][e]   += count * num[e];
						coal_rates_denom[j][e] += count * denom[e];
					}
				}
				if(age_notshared_count[j][bin] > 0){ 
					count = age_notshared_count[j][bin];
					std::fill(num.begin(), num.end(), 0.0);
					std::fill(denom.begin(), denom.end(), 0.0);
					double logl;
					std::vector<double> num2(num_epochs,0.0), denom2(num_epochs,0.0);
					logl = EM[j].EM_notshared_exact(age_bin[bin], num2, denom2);
					logl = EM[j].EM_notshared(age_bin[bin], num, denom);
					assert(logl <= 0.0);
					log_likelihood += count * logl;
					for(int e = 0; e < num_epochs-1; e++){
						//num_diff[e] += count * (num[e]-num2[e]);
						//denom_diff[e] += count * (denom[e]-denom2[e]);
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
					//std::cerr << coal_rates_num[i][e] << " " << coal_rates_denom[i][e] << " " << coal_rates[i][e] << " " << num_diff[e] << " " << denom_diff[e] << std::endl;
				}
			}
		}
		//std::cerr << std::endl;
    //std::cerr << logl_diff << std::endl << std::endl;

		//std::cerr << std::endl;
		os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << "\n";

		//if(log_likelihood - prev_log_likelihood < 0.0) break;

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

