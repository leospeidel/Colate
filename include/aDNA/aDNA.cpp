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
    num_epochs = 5;
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
    //epochs[1] = 1e3/years_per_gen;
    float log_10 = std::log(10);
    for(int e = 1; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;
  }


  ////////////////////////
  //read input sequence (file format? haps/sample? vcf?) and reference sequences (haps/sample? vcf?)
  //For now, assume both are haps/sample and span same positions
  //TODO: make this more flexible

  haps ref_N(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  Data data(ref_N.GetN(), 1);
  ref_N.CloseFile();
  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);

  int i = 0;
  if(options.count("coal") > 0){
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
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    for(int i = 0; i < data.N; i++){
      coal_rates[i].resize(num_epochs);
      coal_rates_num[i].resize(num_epochs);
      coal_rates_denom[i].resize(num_epochs);
      std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
    }
  }

  ////////////////////////////////////////  

  int N = 2;
  int L = mut.info.size();
  int max_iter = 1000;
  int perc = -1;
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
    int snp = 0, snp_input = 0, snp_ref = 0;;
    //iterating over snps in reference
    for(; snp_ref < ref.GetL();){

			//for(int k = 0; k < 1000; k++){
      ref.ReadSNP(sequence_ref, bp_ref);
			snp_ref++;
			//if(snp_ref == ref.GetL()) break;
			//}
			//if(snp_ref == ref.GetL()) break;

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
                    EM[i].EM_shared(mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], num, denom);
										//EM[i].EM_shared(age_begin[snp], age_end[snp], tmrca[snp], num, denom);
                    assert(!std::isnan(denom[0]));
                  }else{
                    EM[i].EM_notshared(mut.info[snp].age_begin, mut.info[snp].age_end, tmrca[snp], num, denom);
										//EM[i].EM_notshared(age_begin[snp], age_end[snp], tmrca[snp], num, denom);
                    assert(!std::isnan(denom[0]));
                  }

                  for(int e = 0; e < num_epochs; e++){
                    coal_rates_num[i][e]   += num[e];
                    coal_rates_denom[i][e] += denom[e];
                  }

                //}
              //}

            }
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

        coal_rates_num[i][e] = 0.0;
        coal_rates_denom[i][e] = 0.0;
      }
    }
    ref.CloseFile();
    input.CloseFile();

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

