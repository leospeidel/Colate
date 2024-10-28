#include "coal.cpp"
//#include "coal_old.cpp"
#include "cxxopts.hpp"
#include <string>

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("CoalRate");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("anc", "Filename of file containing trees.", cxxopts::value<std::string>())
    ("mut", "Filename of file containing mut.", cxxopts::value<std::string>())
		("chr", "Optional: File specifying chromosomes to use.", cxxopts::value<std::string>()) 
    ("bins", "Optional: Epoch boundaries 10^(seq(x,y,stepsize)) [format: x,y,stepsize].", cxxopts::value<std::string>())
    ("years_per_gen", "Optional: Years per generation.", cxxopts::value<float>())
		("seed", "Optional: Seed for random number generator (int)", cxxopts::value<int>())
		("num_bootstraps", "Optional: Number of bootstraps.", cxxopts::value<int>())
    ("poplabels", "Optional: Filename of file containing population labels. One of two formats: Either standard 4 column Relate poplabels format, or local ancestry format. This format lists all possible labels in the first row and in subsequent rows a space-delimited string of chrom BP and integer of the local ancestry label for each haplotype.", cxxopts::value<std::string>())
    ("i,input", "Filename of input.", cxxopts::value<std::string>())
    ("o,output", "Filename of output.", cxxopts::value<std::string>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

	if(!mode.compare("tree")){

		coal(options);

	}else	if(!mode.compare("local_ancestry")){

		coal_localancestry(options);

	}else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "local_ancestry." << std::endl;

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


