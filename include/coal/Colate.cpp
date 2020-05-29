#include "coal.cpp"
#include "coal_old.cpp"
#include "cxxopts.hpp"
#include <string>

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("Colate");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("anc", "Filename of file containing trees.", cxxopts::value<std::string>())
    ("mut", "Filename of file containing mut.", cxxopts::value<std::string>())
		("target_vcf", "Filename of target vcf", cxxopts::value<std::string>())
		("reference_vcf", "Filename of reference vcf", cxxopts::value<std::string>())
    ("target_age", "Target age in generations", cxxopts::value<std::string>())
    ("reference_age", "Reference age in generations", cxxopts::value<std::string>())
		("regularise", "Optional: Specify whether or not to regularise likelihood. 0: no, 1: yes.", cxxopts::value<int>())
		("chr", "Optional: File specifying chromosomes to use.", cxxopts::value<std::string>()) 
    ("bins", "Optional: Epoch boundaries 10^(seq(x,y,stepsize)) [format: x,y,stepsize].", cxxopts::value<std::string>())
    ("coal", "Filename of file containing coalescence rates.", cxxopts::value<std::string>()) 
		("seed", "Optional: Seed for random number generator (int)", cxxopts::value<std::string>())
    ("num_bins", "Optional: Number of bins.", cxxopts::value<int>())
		("correction", "Optional: Specify whether or not to correct for bias due to SNPs segregating. 0: no, 1: yes.", cxxopts::value<int>())
		("haps", "Filename of haps file (Output file format of Shapeit).", cxxopts::value<std::string>())
		("sample", "Filename of sample file (Output file format of Shapeit).", cxxopts::value<std::string>())
    ("i,input", "Filename of input.", cxxopts::value<std::string>())
    ("o,output", "Filename of output.", cxxopts::value<std::string>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

	if(!mode.compare("tree")){

		coal(options);

	}else if(!mode.compare("mut")){

    mut(options);

	}else if(!mode.compare("mut_old")){

    mut_old(options);

	}else if(!mode.compare("mut_tree")){

		mut_tree_fast(options);

	}else if(!mode.compare("mut_perturb")){

		mut_perturb(options);

	}else if(!mode.compare("mut_logl")){

		mut_logl(options);

	}else if(!mode.compare("mut_simple")){

    mut_fast_simplified_all(options);

	}else if(!mode.compare("make_mut_with_out")){

    //make mut
		make_mut_incl_out(options);

	}else if(!mode.compare("test_vcf")){

		//make mut
		test_vcf(options);

	}else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "tree, mut, mut_tree, mut_perturb, mut_logl, mut_simple, make_mut_with_out." << std::endl;

  }

  bool help = false;
  if(!options.count("mode")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

  return 0;

}


