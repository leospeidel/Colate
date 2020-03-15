#include "ftree_func.cpp"
#include "cxxopts.hpp"
#include <string>

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("ftree");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("anc", "Filename of file containing trees.", cxxopts::value<std::string>())
    ("mut", "Filename of file containing mut.", cxxopts::value<std::string>())
    ("haps", "Filename of haps file (Output file format of Shapeit).", cxxopts::value<std::string>())
    ("sample", "Filename of sample file (Output file format of Shapeit).", cxxopts::value<std::string>())
    ("num_bins", "Optional: Number of bins.", cxxopts::value<int>())
    ("coal", "Filename of file containing coalescence rates.", cxxopts::value<std::string>()) 
    ("poplabels", "Filename of poplabels file.", cxxopts::value<std::string>())
    ("i,input", "Filename of input.", cxxopts::value<std::string>())
    ("o,output", "Filename of output.", cxxopts::value<std::string>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("threepop")){

    threepop(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "threepop." << std::endl;

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

