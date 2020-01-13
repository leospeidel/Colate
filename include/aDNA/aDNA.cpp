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


void
MLE_constribution_shared(std::vector<float>& coal_rates, float age_begin, float age_end, float tmrca, std::vector<float>& epochs, std::vector<float>& num, std::vector<float>& denom){

}

void
MLE_constribution_notshared(std::vector<float>& coal_rates, float age_begin, float age_end, float tmrca, std::vector<float>& epochs, std::vector<float>& num, std::vector<float>& denom){

}


void
aDNA(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("poplabels") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, poplabels, output." << std::endl;
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
  std::vector<float> tmrca(mut.info.size());

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

  haps ref(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  Data data(ref.GetN(), ref.GetL());
  std::vector<std::vector<float>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);
  for(int i = 0; i < data.N; i++){
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
  }

  haps input((options["input"].as<std::string>() + ".haps").c_str(), (options["input"].as<std::string>() + ".sample").c_str());
  assert(data.L == input.GetL());

  int bp_ref, bp_input;
  std::vector<char> sequence_ref(data.N), sequence_input(data.N);
  std::vector<float> num(num_epochs), denom(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1
  for(int snp = 0; snp < data.L; snp++){

    ref.ReadSNP(sequence_ref, bp_ref);
    input.ReadSNP(sequence_input, bp_input);
    assert(bp_ref == bp_input);          //TODO: relax this eventually
    assert(bp_ref == mut.info[snp].pos); //TODO: relax this eventually

    //calculate contribution to MLE if shared and non-shared
    for(int i = 0; i < data.N; i++){
      if(sequence_ref[i] == '1'){
      
        if(sequence_input[i] == '0'){
          MLE_constribution_shared(coal_rates[i], mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], epochs, num, denom);
        }else{
          MLE_constribution_notshared(coal_rates[i], mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], epochs, num, denom);
        }
      
        for(int e = 0; e < num_epochs; e++){
          coal_rates_num[i][e]   += num[e];
          coal_rates_denom[i][e] += denom[e];
        }
          
      }
    }

  }

  for(int i = 0; i < data.N; i++){
    for(int e = 0; e < num_epochs; e++){
      coal_rates[i][e] = coal_rates_num[i][e]/coal_rates_denom[i][e];
    }
  }

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

