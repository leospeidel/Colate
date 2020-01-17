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
    std::cout << "Needed: anc, mut, haps, sample, input, output." << std::endl;
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

  //get TMRCA at each SNP
  std::vector<float> tmrca(mut.info.size(), 0.0);
  std::vector<float>::iterator it_tmrca = tmrca.begin();
 
  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
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
        it_tmrca++;
        it_mut++;
      }
    }

  }

  //decide on epochs
  //TODO: custom epoch boundaries
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epochs(num_epochs);
  epochs[0] = 0.0;
  epochs[1] = 1e3/years_per_gen;
  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){
    epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
  }
  epochs[num_epochs-1] = 1e8/years_per_gen;

  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)
  //For now, assume both are haps/sample and span same positions
  //TODO: make this more flexible

  double initial_coal_rate = 1.0/5000.0;

  haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  Data data(ref_N.GetN(), ref_N.GetL());
  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);
  for(int i = 0; i < data.N; i++){
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
  }

  for(int iter = 0; iter < 1; iter++){
  haps ref(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  haps input((options["input"].as<std::string>() + ".haps.gz").c_str(), (options["input"].as<std::string>() + ".sample.gz").c_str());
  assert(data.L == input.GetL());

  int bp_ref, bp_input;
  std::vector<char> sequence_ref(data.N), sequence_input(data.N);
  std::vector<float> num(num_epochs, 0.0), denom(num_epochs, 0.0); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1
  int snp = 0; 
  for(; snp < data.L; snp++){

    ref.ReadSNP(sequence_ref, bp_ref);
    input.ReadSNP(sequence_input, bp_input);
    assert(bp_ref == bp_input);          //TODO: relax this eventually
    while(bp_ref > mut.info[snp].pos){
      snp++;
      if(snp == data.L) break;
    }
    if(snp == data.L) break;
    assert(bp_ref == mut.info[snp].pos); //TODO: relax this eventually

    //calculate contribution to MLE if shared and non-shared
    for(int i = 1; i < data.N; i++){
      if(mut.info[snp].age_end > 0){
        if(sequence_ref[i] == '1'){
        
          if(sequence_input[0] == '1'){
            //std::cerr << "shared " << i << " " << mut.info[snp].age_begin << " " << mut.info[snp].age_end << " " << tmrca[snp] << std::endl;
            EM_shared(coal_rates[i], mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], epochs, num, denom);
          }else{
            //std::cerr << "not shared " << i << " " << mut.info[snp].age_begin << " " << mut.info[snp].age_end << " " << tmrca[snp] << std::endl;
            EM_notshared(coal_rates[i], mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], epochs, num, denom);
          }
        
          for(int e = 0; e < num_epochs; e++){
            coal_rates_num[i][e]   += num[e];
            coal_rates_denom[i][e] += denom[e];
          }
            
        }
      }
    }

  }

  for(int i = 1; i < data.N; i++){
    for(int e = 0; e < num_epochs; e++){
      if(coal_rates_denom[i][e] == 0 && e > 0){
        coal_rates[i][e] = coal_rates[i][e-1];
      }else{
        coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
      }
      //std::cerr << coal_rates[i][e] << " ";
      coal_rates_num[i][e] = 0.0;
      coal_rates_denom[i][e] = 0.0;
    }
    //std::cerr << std::endl;
  }
  //std::cerr << std::endl;

  }

  std::ofstream os(options["output"].as<std::string>());

  for(int i = 0; i < data.N; i++){
    os << i << " ";
  }
  os << std::endl;
 
  for(int e = 0; e < num_epochs; e++){
    os << epochs[e] << " ";
  }
  os << std::endl;

  for(int i = 1; i < data.N; i++){
    os << "0 " << i << " ";
    for(int e = 0; e < num_epochs; e++){
      os << coal_rates[i][e] << " ";
    }
    os << std::endl;
  }

  os.close();

  //should iterate here.

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

