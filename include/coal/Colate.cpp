#include "coal.cpp"
//#include "coal_old.cpp"
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
		("target_bcf", "Filename of target bcf", cxxopts::value<std::string>())
		("reference_bcf", "Filename of reference bcf", cxxopts::value<std::string>())
		("target_mask", "Fasta file containing target mask", cxxopts::value<std::string>())
		("reference_mask", "Fasta file containing reference mask", cxxopts::value<std::string>())
    ("target_table", "Filename of target table containing sampled reads. Columns CHR BP allele", cxxopts::value<std::string>())
		("target_bam", "Filename of target bam", cxxopts::value<std::string>())
		("reference_bam", "Filename of reference bam", cxxopts::value<std::string>())
		("target_tmp", "Filename of target tmp file", cxxopts::value<std::string>())
		("reference_tmp", "Filename of reference tmp file", cxxopts::value<std::string>())
    ("target_age", "Target age in years", cxxopts::value<std::string>())
    ("reference_age", "Reference age in years", cxxopts::value<std::string>())
		("ref_genome", "Reference genome", cxxopts::value<std::string>())
		("anc_genome", "Ancestral genome", cxxopts::value<std::string>())
		("mask", "Genomic mask", cxxopts::value<std::string>())
		("mask_cutoff", "Remove trees with average passing rate below this cutoff. Default 0.9.", cxxopts::value<float>())
		("chr", "Optional: File specifying chromosomes to use.", cxxopts::value<std::string>()) 
    ("bins", "Optional: Epoch boundaries 10^(seq(x,y,stepsize)) [format: x,y,stepsize]. In years.", cxxopts::value<std::string>())
		("lineage_bin", "Optional: Epoch boundary of focal lineage when using CondCoalRates. 10^lineage_bin. In years", cxxopts::value<float>())
		("outgroup_tmrca", "TMRCA to outgroup in years. Default: 10e6 years", cxxopts::value<float>())
    ("years_per_gen", "Optional: Years per generation.", cxxopts::value<float>())
    ("coal", "Filename of file containing coalescence rates.", cxxopts::value<std::string>()) 
		("seed", "Optional: Seed for random number generator (int)", cxxopts::value<int>())
		("num_bootstraps", "Optional: Number of bootstraps.", cxxopts::value<int>())
		("filters", "Optional: Parameters for filtering BAMs. Format: MAPQ,LEN,MAX_MISMATCH_TO_REF. Default: 20,30,10", cxxopts::value<std::string>())
    ("strandfilter", "Optional: When parsing BAMs, filter potential deaminated sites by strand orientation.", cxxopts::value<bool>())
    ("groups", "Names of groups of interest for conditional coalescence rates", cxxopts::value<std::string>()) 
    ("poplabels", "Optional: Filename of file containing population labels. One of two formats: Either standard 4 column Relate poplabels format, or local ancestry format. This format lists all possible labels in the first row and in subsequent rows a space-delimited string of chrom BP and integer of the local ancestry label for each haplotype.", cxxopts::value<std::string>())
		("map", "Genetic map. (Prefix if for multiple chr)", cxxopts::value<std::string>())
    ("i,input", "Filename of input.", cxxopts::value<std::string>())
    ("o,output", "Filename of output.", cxxopts::value<std::string>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("mut")){

    mut(options);
 
	}else if(!mode.compare("preprocess_mut")){

    //make mut
		preprocess_mut(options);

	}else if(!mode.compare("make_tmp")){

    //make mut
	  make_tmp(options);

	}else if(!mode.compare("print_tmp")){

		//make mut
		print_tmp(options);

	}else if(!mode.compare("compare_tmp")){

		//make mut
		compare_tmp(options);

  }else if(!mode.compare("count_topo")){

    //make mut
    count_topo(options);

  }else if(!mode.compare("calc_depth")){

		//make mut
		calc_depth(options);

  }else if(!mode.compare("get_deam")){

    //make mut
    get_deam(options);

  }else if(!mode.compare("CondCoalRates")){

		//make mut
		CondCoalRates(options);

	}else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "preprocess_mut, make_tmp, mut, calc_depth, print_tmp, CondCoalRates." << std::endl;

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


