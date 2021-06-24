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
	if(!options.count("input") || !options.count("output") || !options.count("bins")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output, bins. Optional: years_per_gen, coal" << std::endl;
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

		double min_epoch = 0.0;
		double epoch_lower, epoch_upper, epoch_step;
		double log_10 = log(10);
		double years_per_gen = 28.0;
		if(options.count("years_per_gen")){
			years_per_gen = options["years_per_gen"].as<float>();
		}

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
		epoch_boundary = epoch_lower;
		while(epoch_boundary < epoch_upper){
			epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
			ep++;
			epoch_boundary += epoch_step;
		}
		epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
		epochs.push_back( std::max(1e8, 10*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();	

	}

	int num_bootstrap = 1;
	if(options.count("num_bootstraps")) num_bootstrap = options["num_bootstraps"].as<int>();
	int block_size = 5000;

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

		std::cerr << "CHR " << chromosomes[chr] << ":\n";
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
			ct.populate(mtr.tree, num_bases_tree_persists);	
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
coal_localancestry(cxxopts::Options& options){

	//Program options

	bool help = false;
	if(!options.count("input") || !options.count("output") || !options.count("poplabels") || !options.count("bins")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output, poplabels, bins. Optional: years_per_gen, coal" << std::endl;
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
	std::string line, tmp;
	if(options.count("coal") > 0){
		is.open(options["coal"].as<std::string>());
		getline(is, line);
		getline(is, line);
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
		double epoch_lower, epoch_upper, epoch_step;
		double log_10 = log(10);
		double years_per_gen = 28.0;
		if(options.count("years_per_gen")){
			years_per_gen = options["years_per_gen"].as<float>();
		}

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
		epoch_boundary = epoch_lower;
		while(epoch_boundary < epoch_upper){
			epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
			ep++;
			epoch_boundary += epoch_step;
		}
		epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
		epochs.push_back( std::max(1e8, 10*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();	

	}

	int num_bootstrap = 1;
	if(options.count("num_bootstraps")) num_bootstrap = options["num_bootstraps"].as<int>();
	int block_size = 5000;


	MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut; //iterator for mut file

	std::vector<std::string> chromosomes, filenames;
	if(options.count("chr") > 0){

		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		while(getline(is_chr, line)){
			chromosomes.push_back(line);
			filenames.push_back(options["input"].as<std::string>() + "_chr" + line);
		}
		is_chr.close();

	}else{
		chromosomes.resize(1);
		chromosomes[0] = "NA";
		filenames.push_back(options["input"].as<std::string>());
	}

	//local ancestry
	std::vector<std::string> lchrom, unique_groups;
	std::vector<int> lbp;
	std::vector<std::vector<int>> group;

	igzstream is_assign(options["poplabels"].as<std::string>());
	if(is_assign.fail()){
		std::cerr << "Error: Failed to open " << options["poplabels"].as<std::string>() << std::endl;
		exit(1);
	}
	getline(is_assign, line);
	int i = 0;
	while(i < line.size()){
		tmp.clear();
		while(line[i] != ' ' && line[i] != '\t'){
			tmp += line[i];
			i++;
			if(i == line.size()) break;
		}
		i++;
		unique_groups.push_back(tmp);
	}
	int num_groups = unique_groups.size();

	std::vector<int> group_tmp;
	int val;
	std::string current_chr;
	while(getline(is_assign, line)){

		std::istringstream iss(line);
		if(group_tmp.size() == 0){
			iss >> tmp;
			lchrom.push_back(tmp);
			current_chr = tmp;
			iss >> val;
			lbp.push_back(val);
			if(val != 0){
				std::cerr << "Error: First entry for new chr has to start at BP = 0" << std::endl;
				exit(1);
			}
			while(iss >> val) group_tmp.push_back(val);
		}else{
			iss >> tmp;
			lchrom.push_back(tmp);
			iss >> val;
			lbp.push_back(val);
			if(current_chr != lchrom[lchrom.size()-1] && val != 0){
				std::cerr << "Error: First entry for new chr has to start at BP = 0" << std::endl;
				exit(1);
			}
			current_chr = lchrom[lchrom.size()-1];
			int k = 0;
			while(iss >> group_tmp[k]) k++;
		}
		group.push_back(group_tmp);

	}

	//Calculate coal rates
	coal_LA ct(epochs, num_bootstrap, block_size, num_groups);

	int local_index = 0;
	for(int chr = 0; chr < chromosomes.size(); chr++){

		std::cerr << "CHR " << chromosomes[chr] << ":\n";
		AncMutIterators ancmut(filenames[chr] + ".anc", filenames[chr] + ".mut");
		float num_bases_tree_persists = 0.0;

		if(chromosomes.size() > 1 || chromosomes[chr] != "NA") assert(lchrom[local_index] == chromosomes[chr]);
		assert(lbp[local_index] == 0);

		ct.update_ancmut(ancmut);

		int tree_count = 0, perc = -1;
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		while(num_bases_tree_persists >= 0.0){

			if( (int) (((double)tree_count)/ancmut.NumTrees() * 100.0) > perc ){
				perc = (int) (((double)tree_count)/ancmut.NumTrees() * 100.0);
				std::cerr << "[" << perc << "%]\r";
			}
			tree_count++;

			int bp_start = (*it_mut).pos;
			int t = (*it_mut).tree;
			while((*it_mut).tree == t) it_mut++;
			int bp_end = (*it_mut).pos;
			if(bp_end == bp_start) bp_end++;

			assert(group[local_index].size() == ancmut.NumTips());
			if(local_index < group.size()-1){

				if(bp_end > lbp[local_index+1] && (lchrom[local_index+1] == chromosomes[chr] || chromosomes[chr] == "NA") ){

					double frac = (lbp[local_index+1] - bp_start)/(bp_end-bp_start);
					assert(frac <= 1.0);
					assert(frac >= 0.0);
					ct.populate(mtr.tree, num_bases_tree_persists * frac, group[local_index], true);
					local_index++;
					if(local_index == group.size()){
						local_index--;
					}
					while(bp_end > lbp[local_index+1] && (lchrom[local_index+1] == chromosomes[chr] || chromosomes[chr] == "NA") ){
						frac = (lbp[local_index+1] - lbp[local_index])/(bp_end-bp_start);
						assert(frac <= 1.0);
						assert(frac >= 0.0);
						ct.populate(mtr.tree, num_bases_tree_persists * frac, group[local_index], false);
						local_index++;
						if(local_index == group.size()){
							local_index--;
							break;
						}
					}
					frac = (bp_end - lbp[local_index])/(bp_end-bp_start);
					assert(frac <= 1.0);
					assert(frac >= 0.0);
					ct.populate(mtr.tree, num_bases_tree_persists * frac, group[local_index], false);
				}else{
					ct.populate(mtr.tree, num_bases_tree_persists, group[local_index], true);
				}

			}else{
				ct.populate(mtr.tree, num_bases_tree_persists, group[local_index], true);
			}

			num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		}
		local_index++;

		std::cerr << std::endl;

	}

	ct.Dump(options["output"].as<std::string>() + ".coal", unique_groups);

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


//////////////////////
int
parse_vcf(std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_target_mask, std::vector<std::string>& filename_ref_genome, double age, double ref_age, double C, std::mt19937& rng, int num_bases_per_block, std::vector<std::vector<double>>& age_shared_count, std::vector<std::vector<double>>& age_notshared_count, std::vector<std::vector<double>>& age_shared_emp, std::vector<std::vector<double>>& age_notshared_emp){

	std::uniform_real_distribution<double> dist_unif(0,1);
	int num_age_bins = ((int) (log(1e8) * C))+1;
	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0;
	bool target_flip = false;
	std::string ancestral, derived, ancestral_target, derived_target;
	Muts::iterator it_mut;

	float num_samples = 100;
	int DAF_target, DAF_ref;
	int bin_index;

	bool has_mask = false;
	if(filename_target_mask.size() > 0) has_mask = true;
	bool has_ref_genome = false; //new
	if(filename_ref_genome.size() > 0) has_ref_genome = true; //new

	double proportion_with_age = 0.0, total_mut = 0.0;

	std::vector<std::vector<double>>::iterator it_age_shared_count    = age_shared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_shared_emp      = age_shared_emp.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_emp   = age_notshared_emp.begin();
	int current_block_base = 0;
	int num_blocks = 0;

	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

		Mutations mut;
		mut.Read(filename_mut[chr]);
		vcf_parser target(filename_target[chr]);
		fasta mask;
		if(has_mask) mask.Read(filename_target_mask[chr]);
		fasta ref_genome; //new
		if(has_ref_genome) ref_genome.Read(filename_ref_genome[chr]); //new

		current_block_base = 0;
		bp_target = 0, bp_mut = 0;

		target.read_snp();
		bp_target = target.rec->pos + 1;
		target.extract_GT();
		N_target    = target.n * target.ploidy;
		total_mut++;

		double DAF_check = 0;
		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

				if((*it_mut).age_end >= age || (*it_mut).age_end >= age){

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

						bool use = true;          
						if(has_mask && bp_mut < mask.seq.size()){
							if(mask.seq[bp_mut-1] != 'P') use = false;
						}
						if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
						if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;

						if(use){

							//check if mutation exists in ref
							if(bp_target < bp_mut){
								while(target.read_snp() == 0){
									bp_target = target.rec->pos + 1;
									target.extract_GT(); 
									DAF_check = 0;
									for(int i = 0; i < target.n*target.ploidy; i++){
										DAF_check += bcf_gt_allele(target.gt[i]);
									}
									if(DAF_check > 0) total_mut++;
									if(bp_target >= bp_mut) break;
								}
							}
							if(bp_target == bp_mut){
								if(DAF_check > 0) proportion_with_age++;
							}

							DAF_target = 0;
							DAF_ref    = 0;
							if(bp_target == bp_mut){ //exists

								ancestral_target = target.rec->d.allele[0];
								derived_target   = target.rec->d.allele[1];

								bool accept = false, fixed_for_ref = false;
								target_flip = false;
								if(ancestral_target == derived && derived_target == ""){

									//possibility of sharing
									target.extract_GT();

									if(target.n*target.ploidy < 2){
										std::cerr << "Need at least one diploid genome." << std::endl;
										exit(1);
									}

									N_ref = 0;
									N_target = 0;

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
										if(DAF_target != 0 || DAF_ref != 0){
											use = false;
										}else{
											DAF_target = N_target;
											DAF_ref    = N_ref;
										}
									}else{
										use = false;
									}

								}else if( (ancestral_target == ancestral && derived_target == derived) || (ancestral_target == derived && derived_target == ancestral) ){

									target_flip = false;
									if( (ancestral_target == derived && derived_target == ancestral) ) target_flip = true;

									double age_begin = (*it_mut).age_begin;

									//possibility of sharing
									target.extract_GT();

									if(target.n*target.ploidy < 2){
										std::cerr << "Need at least one diploid genome." << std::endl;
										exit(1);
									}

									N_ref = 0;
									N_target = 0;

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
									}else{
										use = false;
									}
								}else{
									use = false;
								}
							}else{ //doesn't exist so use ref genome

								if(has_ref_genome){
									derived_target = ref_genome.seq[bp_mut-1];
									if(derived == derived_target){
										//std::cerr << bp_mut << " " << ancestral << " " << derived << " " << derived_target << " " << (*it_mut).age_begin << " " << (*it_mut).age_end << std::endl;
										DAF_ref = N_ref;
										DAF_target = N_target;
									}else{
										use = false;
									}
								}else{
									use = false;
								}

							}

							if(DAF_ref == 0) use = false;

							if(use){

								while(current_block_base + num_bases_per_block < bp_mut){
									current_block_base += num_bases_per_block;
									it_age_shared_emp++;
									it_age_notshared_emp++;
									it_age_shared_count++;
									it_age_notshared_count++;
									num_blocks++;
								}

								double age_begin = (*it_mut).age_begin;
								if(age_begin < ref_age) age_begin = ref_age;
								num_used_snps++;

								bool skip = false;
								if(age_begin <= age){

									double age_begin2 = age_begin;
									age_begin2 = 0.0;

									int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
									int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
									bin_index  = bin_index1 * num_age_bins + bin_index2;
									(*it_age_shared_emp)[bin_index]    += DAF_target * DAF_ref/((double)N_ref);
									(*it_age_notshared_emp)[bin_index] += (N_target - DAF_target) * DAF_ref/((double)N_ref);

									if(1){
										int j = 0;
										while(j < num_samples){
											skip = false;
											double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
											//if(sampled_age < age) skip = true;
											if(sampled_age < age) sampled_age = age;
											int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
											bin_index = bin_index_age;// * num_age_bins + bin_index_age;

											if(!skip){
												//age_shared_count[bin_index]    += DAF_target * DAF_ref/((double) num_samples);
												(*it_age_notshared_count)[bin_index] += (N_target - DAF_target) * DAF_ref/((double) num_samples * N_ref);
											}
											j++;
										}
									}

								}else{

									int bin_index1  = std::max(0, (int)std::round(log(10*(*it_mut).age_begin)*C)+1);
									double age_end2 = (*it_mut).age_end;
									int bin_index2  = std::max(0, (int)std::round(log(10*age_end2)*C)+1);
									bin_index       = bin_index1 * num_age_bins + bin_index2;

									int j = 0;
									while(j < num_samples){
										skip = false;
										double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
										if(sampled_age < age) skip = true;
										int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
										bin_index = bin_index_age;// * num_age_bins + bin_index_age;

										if(!skip){
											(*it_age_shared_count)[bin_index]    += DAF_target * DAF_ref/((double) num_samples * N_ref);
											(*it_age_notshared_count)[bin_index] += (N_target - DAF_target) * DAF_ref/((double) num_samples * N_ref);
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
		it_age_shared_count++;
		it_age_notshared_count++;
		it_age_shared_emp++;
		it_age_notshared_emp++;
		num_blocks++;
	}

	//std::cerr << proportion_with_age/total_mut << std::endl;
	//std::cerr << num_used_snps << std::endl << std::endl;

	age_shared_count.resize(num_blocks);
	age_notshared_count.resize(num_blocks);
	age_shared_emp.resize(num_blocks);
	age_notshared_emp.resize(num_blocks);

	return num_blocks;

}

int
parse_vcfvcf(std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_ref, std::vector<std::string>& filename_target_mask, std::vector<std::string>& filename_reference_mask, std::vector<std::string>& filename_ref_genome, double age, double ref_age, double C, std::mt19937& rng, int num_bases_per_block, std::vector<std::vector<double>>& age_shared_count, std::vector<std::vector<double>>& age_notshared_count, std::vector<std::vector<double>>& age_shared_emp, std::vector<std::vector<double>>& age_notshared_emp){

	std::uniform_real_distribution<double> dist_unif(0,1);
	int num_age_bins = ((int) (log(1e8) * C))+1;
	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0;
	bool ref_flip = false, target_flip = false;
	std::string ancestral, derived, ancestral_ref, derived_ref, ancestral_target, derived_target;
	Muts::iterator it_mut;

	float num_samples = 100;
	int DAF_target, DAF_ref;
	int bin_index;

	bool has_tar_mask = false, has_ref_mask = false;
	if(filename_target_mask.size() > 0) has_tar_mask = true;
	if(filename_reference_mask.size() > 0) has_ref_mask = true;
	bool has_ref_genome = false; //new
	if(filename_ref_genome.size() > 0) has_ref_genome = true; //new

	std::vector<std::vector<double>>::iterator it_age_shared_count    = age_shared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_shared_emp      = age_shared_emp.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_emp   = age_notshared_emp.begin();
	int current_block_base = 0;
	int num_blocks = 0;

	bool shared = false;
	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

		Mutations mut;
		mut.Read(filename_mut[chr]);
		vcf_parser target(filename_target[chr]);
		vcf_parser ref(filename_ref[chr]);
		fasta tar_mask, ref_mask;
		if(has_tar_mask) tar_mask.Read(filename_target_mask[chr]);
		if(has_ref_mask) ref_mask.Read(filename_reference_mask[chr]);
		fasta ref_genome; //new
		if(has_ref_genome) ref_genome.Read(filename_ref_genome[chr]); //new
		current_block_base = 0;
		bp_ref = 0, bp_target = 0, bp_mut = 0;

		ref.read_snp();
		bp_ref = ref.rec->pos + 1;
		ref.extract_GT();
		N_ref       = ref.n * ref.ploidy; 

		target.read_snp();
		bp_target = target.rec->pos + 1;
		target.extract_GT();
		N_target    = target.n * target.ploidy;

		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			shared = false;
			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

				if((*it_mut).age_end >= age || (*it_mut).age_end >= age){

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

						bool use = true;
						if(has_tar_mask && bp_mut < tar_mask.seq.size()){
							if(tar_mask.seq[bp_mut-1] != 'P') use = false;
						}
						if(has_ref_mask && bp_mut < ref_mask.seq.size()){
							if(ref_mask.seq[bp_mut-1] != 'P') use = false;
						}
						if(ancestral == derived) use = false;
						if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
						if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;

						if(use){
							//I can use this mutation subject to data in target and ref

							//check if mutation exists in ref
							if(bp_ref < bp_mut){
								while(ref.read_snp() == 0){
									bp_ref = ref.rec->pos + 1;
									if(bp_ref >= bp_mut) break;
								}
							}

							DAF_ref = 0;
							if(bp_ref == bp_mut){ //exists

								ancestral_ref = ref.rec->d.allele[0];
								derived_ref   = ref.rec->d.allele[1];
								if( (ancestral_ref == ancestral && derived_ref == derived) || (ancestral_ref == derived && derived_ref == ancestral) ){

									ref_flip = false;
									if( (ancestral_ref == derived && derived_ref == ancestral) ) ref_flip = true;

									ref.extract_GT();
									N_ref       = ref.n * ref.ploidy; 
									bool biallelic = true;
									for(int i = 0; i < ref.n; i++){
										for(int j = 0; j < ref.ploidy; j++){
											DAF_ref += bcf_gt_allele(ref.gt[ref.ploidy*i+j]);
											if(bcf_gt_allele(ref.gt[ref.ploidy*i+j]) > 1) biallelic = false;
										}
									}

									if(biallelic){
										assert(DAF_ref <= N_ref);
										if(ref_flip){
											DAF_ref = N_ref - DAF_ref;
										}
									}else{
										use = false;
									}

								}else{
									use = false;
								}

							}else{ //doesn't exist

								if(has_ref_genome){

									derived_ref = ref_genome.seq[bp_mut-1];
									if(derived == derived_ref){
										DAF_ref = N_ref;
									}else{
										use = false;
									}
								}else{
									use = false;
								}

							}
							if(DAF_ref == 0) use = false;

							if(use){

								//check if mutation exists in target
								if(bp_target < bp_mut){
									while(target.read_snp() == 0){
										bp_target = target.rec->pos + 1;
										if(bp_target >= bp_mut) break;
									}
								} 

								DAF_target = 0;
								if(bp_target == bp_mut){ //exists

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

									if(!accept) use = false;
									if(target_flip){
										DAF_target = N_target - DAF_target;
									}

								}else{

									if(has_ref_genome){
										derived_target = ref_genome.seq[bp_mut-1];
										if(derived == derived_target){
											DAF_target = N_target;
											shared = true;
										}else if(ancestral == derived_target){
											DAF_target = 0;
										}else{
											use = false;
										}
									}else{
										use = false;
									}

								}

							}

						}

						if(use){

							while(current_block_base + num_bases_per_block < bp_mut){
								current_block_base += num_bases_per_block;
								it_age_shared_emp++;
								it_age_notshared_emp++;
								it_age_shared_count++;
								it_age_notshared_count++;
								num_blocks++;
							}

							double age_begin = (*it_mut).age_begin;
							if(age_begin < ref_age) age_begin = ref_age;
							num_used_snps++;

							bool skip = false;
							if(age_begin <= age){

								double age_begin2 = age_begin;
								age_begin2 = 0.0;
								int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
								int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
								bin_index  = bin_index1 * num_age_bins + bin_index2;
								(*it_age_shared_emp)[bin_index]    += DAF_target * DAF_ref/((double) N_ref);
								(*it_age_notshared_emp)[bin_index] += (N_target - DAF_target) * DAF_ref/((double) N_ref);

								int j = 0;
								//if(DAF_target > 0) age_begin = (*it_mut).age_end;
								while(j < num_samples){
									skip = false;
									double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
									//if(sampled_age < age) skip = true;
									if(sampled_age < age) sampled_age = age;
									int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
									bin_index = bin_index_age;// * num_age_bins + bin_index_age;

									if(!skip){
										//age_shared_count[bin_index]    += DAF_target * DAF_ref/num_samples;
										(*it_age_notshared_count)[bin_index] += (N_target - DAF_target) * DAF_ref/((double) N_ref * num_samples);
										//j++;
									}
									j++;
								}

							}else{

								int j = 0;
								while(j < num_samples){
									skip = false;
									double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
									if(sampled_age < age) skip = true;
									int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
									bin_index = bin_index_age;// * num_age_bins + bin_index_age;

									if(!skip){
										(*it_age_shared_count)[bin_index]    += DAF_target * DAF_ref/((double) N_ref * num_samples);
										(*it_age_notshared_count)[bin_index] += (N_target - DAF_target) * DAF_ref/((double) N_ref * num_samples);
										j++;
									}
									//j++;
								}

							} 

						}
					}
				}
			}
		}
		it_age_shared_count++;
		it_age_notshared_count++;
		it_age_shared_emp++;
		it_age_notshared_emp++;
		num_blocks++;

	}

	std::cerr << num_used_snps << std::endl << std::endl;

	age_shared_count.resize(num_blocks);
	age_notshared_count.resize(num_blocks);
	age_shared_emp.resize(num_blocks);
	age_notshared_emp.resize(num_blocks);

	return num_blocks;

}

int
parse_bamvcf(std::string& params, std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_ref, std::vector<std::string>& filename_target_mask, std::vector<std::string>& filename_reference_mask, std::vector<std::string>& filename_ref_genome, double age, double ref_age, double C, std::mt19937& rng, int num_bases_per_block, std::vector<std::vector<double>>& age_shared_count, std::vector<std::vector<double>>& age_notshared_count, std::vector<std::vector<double>>& age_shared_emp, std::vector<std::vector<double>>& age_notshared_emp){

	std::uniform_real_distribution<double> dist_unif(0,1);
	int num_age_bins = ((int) (log(1e8) * C))+1;
	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0, num_used_snps2 = 0;
	bool ref_flip = false, target_flip = false;
	std::string ancestral, derived, ancestral_ref, derived_ref, ancestral_target, derived_target;
	Muts::iterator it_mut;

	float num_samples = 100;
	int DAF_target, AAF_target, DAF_ref;
	int bin_index;

	bool has_tar_mask = false, has_ref_mask = false;
	if(filename_target_mask.size() > 0) has_tar_mask = true;
	if(filename_reference_mask.size() > 0) has_ref_mask = true;
	bool has_ref_genome = false; //new
	if(filename_ref_genome.size() > 0) has_ref_genome = true; //new

	std::vector<std::vector<double>>::iterator it_age_shared_count    = age_shared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_shared_emp      = age_shared_emp.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_emp   = age_notshared_emp.begin();
	int current_block_base = 0;
	int num_blocks = 0;

	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

		Mutations mut;
		mut.Read(filename_mut[chr]);
		bam_parser target(filename_target[chr], params, filename_ref_genome[chr]);
		vcf_parser ref(filename_ref[chr]);	

		fasta tar_mask, ref_mask;
		if(has_tar_mask) tar_mask.Read(filename_target_mask[chr]);
		if(has_ref_mask) ref_mask.Read(filename_reference_mask[chr]);
		fasta ref_genome; //new
		if(has_ref_genome) ref_genome.Read(filename_ref_genome[chr]); //new

		current_block_base = 0;
		bp_ref = 0, bp_target = 0, bp_mut = 0;
		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

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

					bool use = true;
					if(has_tar_mask && bp_mut < tar_mask.seq.size()){
						if(tar_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(has_ref_mask && bp_mut < ref_mask.seq.size()){
						if(ref_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;


					if(use){

						DAF_ref    = 0;
						DAF_target = 0;
						AAF_target = 0;

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


								ref.extract_GT();
								N_ref       = ref.n * ref.ploidy; 
								bool biallelic = true;
								for(int i = 0; i < ref.n; i++){
									for(int j = 0; j < ref.ploidy; j++){
										DAF_ref += bcf_gt_allele(ref.gt[ref.ploidy*i+j]);
										if(bcf_gt_allele(ref.gt[ref.ploidy*i+j]) > 1) biallelic = false;
									}
								}

								if(biallelic){
									assert(DAF_ref <= N_ref);
									if(ref_flip){
										DAF_ref = N_ref - DAF_ref;
									}
								}else{
									use = false;
								}

							}else{
								use = false;
							}

						}else{ //doesn't exist

							if(1){
								if(has_ref_genome){
									derived_ref   = ref_genome.seq[bp_mut-1];
									if(derived == derived_ref){
										DAF_ref = N_ref;
									}else{
										use = false;
									}
								}else{
									use = false;
								}
							}

						}
						if(DAF_ref == 0) use = false;

						if(use){

							//check if mutation exists in target     
							int bp_target = bp_mut - 1;
							target.read_to_pos(bp_target);
							int num_alleles = 0;
							if(target.pos_of_entry[bp_target % target.num_entries] == bp_target){
								int num_reads = 0;
								for(int i = 0; i < 4; i++){
									num_reads += target.count_alleles[bp_target % target.num_entries][i];
									num_alleles += (target.count_alleles[bp_target % target.num_entries][i] > 0);
								}

								if(num_reads > 0){

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
										if(!(num_alleles <= 2)) use = false;
									}else{
										use = false;
									}

								}else{
									use = false;
								}

							}else{
								use = false;
							}
						} 

						//populate age_shared_count and age_notshared_count
						if(use){

							double age_begin = (*it_mut).age_begin;
							if(age_begin < ref_age) age_begin = ref_age;

							while(current_block_base + num_bases_per_block < bp_mut){
								current_block_base += num_bases_per_block;
								it_age_shared_emp++;
								it_age_notshared_emp++;
								it_age_shared_count++;
								it_age_notshared_count++;
								num_blocks++;
							}

							bool skip = false;
							if(age_begin <= age){

								if(1){
									double age_begin2 = age_begin;
									age_begin2 = 0.0;
									int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
									int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
									bin_index  = bin_index1 * num_age_bins + bin_index2;
									(*it_age_shared_emp)[bin_index]    += DAF_target * DAF_ref/((double)N_ref);
									(*it_age_notshared_emp)[bin_index] += AAF_target * DAF_ref/((double)N_ref);

									if(1){
										int j = 0;
										while(j < num_samples){
											skip = false;
											double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
											//if(sampled_age < age) skip = true;
											if(sampled_age < age) sampled_age = age;
											int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
											bin_index = bin_index_age;// * num_age_bins + bin_index_age;

											if(!skip){
												(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
												//j++;
											}
											j++;
										}
									}
								}

							}else{

								int j = 0;
								while(j < num_samples){
									skip = false;
									double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
									if(sampled_age < age) skip = true;
									int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
									bin_index = bin_index_age;// * num_age_bins + bin_index_age;

									if(!skip){
										(*it_age_shared_count)[bin_index]    += DAF_target * DAF_ref/((double)N_ref * num_samples);
										(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
										j++;
									}
								}

							}


						}

					}

				}
			}
		}

		it_age_shared_count++;
		it_age_notshared_count++;
		it_age_shared_emp++;
		it_age_notshared_emp++;
		num_blocks++;
		std::cerr << "Coverage: " << ((double) target.coverage)/target.ref_genome.seq.size() << " " << ((double) target.coverage_after_filter)/target.ref_genome.seq.size() << std::endl; 
	}

	age_shared_count.resize(num_blocks);
	age_notshared_count.resize(num_blocks);
	age_shared_emp.resize(num_blocks);
	age_notshared_emp.resize(num_blocks);

	return num_blocks;

}

int
parse_onebamvcf(std::string& params, std::vector<std::string>& filename_chr, std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_ref, std::vector<std::string>& filename_target_mask, std::vector<std::string>& filename_reference_mask, std::vector<std::string>& filename_ref_genome, double age, double ref_age, double C, std::mt19937& rng, int num_bases_per_block, std::vector<std::vector<double>>& age_shared_count, std::vector<std::vector<double>>& age_notshared_count, std::vector<std::vector<double>>& age_shared_emp, std::vector<std::vector<double>>& age_notshared_emp){

	std::uniform_real_distribution<double> dist_unif(0,1);
	int num_age_bins = ((int) (log(1e8) * C))+1;
	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0, num_used_snps2 = 0;
	bool ref_flip = false, target_flip = false;
	std::string ancestral, derived, ancestral_ref, derived_ref, ancestral_target, derived_target;
	Muts::iterator it_mut;

	float num_samples = 100;
	int DAF_target, AAF_target, DAF_ref;
	int bin_index;

	bool has_tar_mask = false, has_ref_mask = false;
	if(filename_target_mask.size() > 0) has_tar_mask = true;
	if(filename_reference_mask.size() > 0) has_ref_mask = true;
	bool has_ref_genome = false; //new
	if(filename_ref_genome.size() > 0) has_ref_genome = true; //new

	assert(filename_target.size() == 1);
	bam_parser target(filename_target[0], params);

	std::vector<std::vector<double>>::iterator it_age_shared_count    = age_shared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_shared_emp      = age_shared_emp.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_emp   = age_notshared_emp.begin();
	int current_block_base = 0;
	int num_blocks = 0;

	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

		target.assign_contig(filename_chr[chr], filename_ref_genome[chr]);

		Mutations mut;
		mut.Read(filename_mut[chr]);
		vcf_parser ref(filename_ref[chr]);	

		fasta tar_mask, ref_mask;
		if(has_tar_mask) tar_mask.Read(filename_target_mask[chr]);
		if(has_ref_mask) ref_mask.Read(filename_reference_mask[chr]);
		fasta ref_genome; //new
		if(has_ref_genome) ref_genome.Read(filename_ref_genome[chr]); //new

		current_block_base = 0;
		bp_ref = 0, bp_target = 0, bp_mut = 0;
		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

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

					bool use = true;
					if(has_tar_mask && bp_mut < tar_mask.seq.size()){
						if(tar_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(has_ref_mask && bp_mut < ref_mask.seq.size()){
						if(ref_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;


					if(use){

						DAF_ref    = 0;
						DAF_target = 0;
						AAF_target = 0;

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


								ref.extract_GT();
								N_ref       = ref.n * ref.ploidy; 
								bool biallelic = true;
								for(int i = 0; i < ref.n; i++){
									for(int j = 0; j < ref.ploidy; j++){
										DAF_ref += bcf_gt_allele(ref.gt[ref.ploidy*i+j]);
										if(bcf_gt_allele(ref.gt[ref.ploidy*i+j]) > 1) biallelic = false;
									}
								}

								if(biallelic){
									assert(DAF_ref <= N_ref);
									if(ref_flip){
										DAF_ref = N_ref - DAF_ref;
									}
								}else{
									use = false;
								}

							}else{
								use = false;
							}

						}else{ //doesn't exist

							if(1){
								if(has_ref_genome){
									derived_ref   = ref_genome.seq[bp_mut-1];
									if(derived == derived_ref){
										DAF_ref = N_ref;
									}else{
										use = false;
									}
								}else{
									use = false;
								}
							}

						}
						if(DAF_ref == 0) use = false;

						if(use){

							//check if mutation exists in target     
							int bp_target = bp_mut - 1;
							target.read_to_pos(bp_target);
							int num_alleles = 0;
							if(target.pos_of_entry[bp_target % target.num_entries] == bp_target){
								int num_reads = 0;
								for(int i = 0; i < 4; i++){
									num_reads += target.count_alleles[bp_target % target.num_entries][i];
									num_alleles += (target.count_alleles[bp_target % target.num_entries][i] > 0);
								}

								if(num_reads > 0){

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
										if(!(num_alleles <= 2)) use = false;
									}else{
										use = false;
									}

								}else{
									use = false;
								}

							}else{
								use = false;
							}
						} 

						//populate age_shared_count and age_notshared_count
						if(use){

							double age_begin = (*it_mut).age_begin;
							if(age_begin < ref_age) age_begin = ref_age;
							while(current_block_base + num_bases_per_block < bp_mut){
								current_block_base += num_bases_per_block;
								it_age_shared_emp++;
								it_age_notshared_emp++;
								it_age_shared_count++;
								it_age_notshared_count++;
								num_blocks++;
							}

							bool skip = false;
							if(age_begin <= age){

								if(1){
									double age_begin2 = age_begin;
									age_begin2 = 0.0;
									int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
									int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
									bin_index  = bin_index1 * num_age_bins + bin_index2;
									(*it_age_shared_emp)[bin_index]    += DAF_target * DAF_ref/((double)N_ref);
									(*it_age_notshared_emp)[bin_index] += AAF_target * DAF_ref/((double)N_ref);

									if(1){
										int j = 0;
										while(j < num_samples){
											skip = false;
											double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
											//if(sampled_age < age) skip = true;
											if(sampled_age < age) sampled_age = age;
											int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
											bin_index = bin_index_age;// * num_age_bins + bin_index_age;

											if(!skip){
												(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
												//j++;
											}
											j++;
										}
									}
								}

							}else{

								int j = 0;
								while(j < num_samples){
									skip = false;
									double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
									if(sampled_age < age) skip = true;
									int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
									bin_index = bin_index_age;// * num_age_bins + bin_index_age;

									if(!skip){
										(*it_age_shared_count)[bin_index]    += DAF_target * DAF_ref/((double)N_ref * num_samples);
										(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
										j++;
									}
								}

							}


						}

					}

				}
			}
		}

		it_age_shared_count++;
		it_age_notshared_count++;
		it_age_shared_emp++;
		it_age_notshared_emp++;
		num_blocks++;
		std::cerr << "Coverage: " << ((double) target.coverage)/target.ref_genome.seq.size() << " " << ((double) target.coverage_after_filter)/target.ref_genome.seq.size() << std::endl; 
	}

	age_shared_count.resize(num_blocks);
	age_notshared_count.resize(num_blocks);
	age_shared_emp.resize(num_blocks);
	age_notshared_emp.resize(num_blocks);

	return num_blocks;

}

int
parse_onebambam(std::string& params, std::vector<std::string>& filename_chr, std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_ref, std::vector<std::string>& filename_target_mask, std::vector<std::string>& filename_reference_mask, std::vector<std::string>& filename_ref_genome, double age, double ref_age, double C, std::mt19937& rng, int num_bases_per_block, std::vector<std::vector<double>>& age_shared_count, std::vector<std::vector<double>>& age_notshared_count, std::vector<std::vector<double>>& age_shared_emp, std::vector<std::vector<double>>& age_notshared_emp){

	std::uniform_real_distribution<double> dist_unif(0,1);
	int num_age_bins = ((int) (log(1e8) * C))+1;
	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0, num_used_snps2 = 0;
	bool ref_flip = false, target_flip = false;
	std::string ancestral, derived, ancestral_ref, derived_ref, ancestral_target, derived_target;
	Muts::iterator it_mut;

	float num_samples = 100;
	int DAF_target, AAF_target, DAF_ref, AAF_ref;
	int bin_index;

	bool has_tar_mask = false, has_ref_mask = false;
	if(filename_target_mask.size() > 0) has_tar_mask = true;
	if(filename_reference_mask.size() > 0) has_ref_mask = true;

	assert(filename_target.size() == 1);
	bam_parser target(filename_target[0], params);
	assert(filename_ref.size() == 1);
	bam_parser reference(filename_ref[0], params);

	std::vector<std::vector<double>>::iterator it_age_shared_count    = age_shared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_shared_emp      = age_shared_emp.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_emp   = age_notshared_emp.begin();
	int current_block_base = 0;
	int num_blocks = 0;

	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

		target.assign_contig(filename_chr[chr], filename_ref_genome[chr]);
		reference.assign_contig(filename_chr[chr], filename_ref_genome[chr]);

		Mutations mut;
		mut.Read(filename_mut[chr]);

		fasta tar_mask, ref_mask;
		if(has_tar_mask) tar_mask.Read(filename_target_mask[chr]);
		if(has_ref_mask) ref_mask.Read(filename_reference_mask[chr]);

		current_block_base = 0;
		bp_ref = 0, bp_target = 0, bp_mut = 0;
		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

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

					bool use = true;
					if(has_tar_mask && bp_mut < tar_mask.seq.size()){
						if(tar_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(has_ref_mask && bp_mut < ref_mask.seq.size()){
						if(ref_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;

					if(use){

						DAF_ref    = 0;
						AAF_ref    = 0;
						DAF_target = 0;
						AAF_target = 0;

						//check if bp exists in ref
						int bp_ref = bp_mut - 1;
						reference.read_to_pos(bp_ref);
						int num_alleles = 0;
						if(reference.pos_of_entry[bp_ref % reference.num_entries] == bp_ref){
							int num_reads = 0;
							for(int i = 0; i < 4; i++){
								num_reads += reference.count_alleles[bp_ref % reference.num_entries][i];
								num_alleles += (reference.count_alleles[bp_ref % reference.num_entries][i] > 0);
							}

							if(num_reads > 0){

								if(ancestral == "A"){
									AAF_ref = reference.count_alleles[bp_ref % reference.num_entries][0];
								}else if(ancestral == "C"){
									AAF_ref = reference.count_alleles[bp_ref % reference.num_entries][1];
								}else if(ancestral == "G"){
									AAF_ref = reference.count_alleles[bp_ref % reference.num_entries][2];
								}else if(ancestral == "T"){
									AAF_ref = reference.count_alleles[bp_ref % reference.num_entries][3];
								}
								if(derived == "A"){
									DAF_ref = reference.count_alleles[bp_ref % reference.num_entries][0];
								}else if(derived == "C"){
									DAF_ref = reference.count_alleles[bp_ref % reference.num_entries][1];
								}else if(derived == "G"){
									DAF_ref = reference.count_alleles[bp_ref % reference.num_entries][2];
								}else if(derived == "T"){
									DAF_ref = reference.count_alleles[bp_ref % reference.num_entries][3];
								}

								if( (AAF_ref > 0 || DAF_ref > 0)){
									if(!(num_alleles <= 2)) use = false;
								}else{
									use = false;
								}

							}else{
								use = false;
							}

						}else{
							use = false;
						}
					} 
					if(DAF_ref == 0) use = false;
					N_ref = DAF_ref + AAF_ref;

					if(use){

						//check if mutation exists in target     
						int bp_target = bp_mut - 1;
						target.read_to_pos(bp_target);
						int num_alleles = 0;
						if(target.pos_of_entry[bp_target % target.num_entries] == bp_target){
							int num_reads = 0;
							for(int i = 0; i < 4; i++){
								num_reads += target.count_alleles[bp_target % target.num_entries][i];
								num_alleles += (target.count_alleles[bp_target % target.num_entries][i] > 0);
							}

							if(num_reads > 0){

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
									if(!(num_alleles <= 2)) use = false;
								}else{
									use = false;
								}

							}else{
								use = false;
							}

						}else{
							use = false;
						}
					} 

					//populate age_shared_count and age_notshared_count
					if(use){

						double age_begin = (*it_mut).age_begin;
						if(age_begin < ref_age) age_begin = ref_age;
						while(current_block_base + num_bases_per_block < bp_mut){
							current_block_base += num_bases_per_block;
							it_age_shared_emp++;
							it_age_notshared_emp++;
							it_age_shared_count++;
							it_age_notshared_count++;
							num_blocks++;
						}

						bool skip = false;
						if(age_begin <= age){

							if(1){
								double age_begin2 = age_begin;
								age_begin2 = 0.0;
								int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
								int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
								bin_index  = bin_index1 * num_age_bins + bin_index2;
								(*it_age_shared_emp)[bin_index]    += DAF_target * DAF_ref/((double)N_ref);
								(*it_age_notshared_emp)[bin_index] += AAF_target * DAF_ref/((double)N_ref);

								if(1){
									int j = 0;
									while(j < num_samples){
										skip = false;
										double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
										//if(sampled_age < age) skip = true;
										if(sampled_age < age) sampled_age = age;
										int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
										bin_index = bin_index_age;// * num_age_bins + bin_index_age;

										if(!skip){
											(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
											//j++;
										}
										j++;
									}
								}
							}

						}else{

							int j = 0;
							while(j < num_samples){
								skip = false;
								double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
								if(sampled_age < age) skip = true;
								int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
								bin_index = bin_index_age;// * num_age_bins + bin_index_age;

								if(!skip){
									(*it_age_shared_count)[bin_index]    += DAF_target * DAF_ref/((double)N_ref * num_samples);
									(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
									j++;
								}
							}

						}


					}

				}
			}
		}

		it_age_shared_count++;
		it_age_notshared_count++;
		it_age_shared_emp++;
		it_age_notshared_emp++;
		num_blocks++;
		std::cerr << "Target coverage: " << ((double) target.coverage)/target.ref_genome.seq.size() << " " << ((double) target.coverage_after_filter)/target.ref_genome.seq.size() << std::endl;
		std::cerr << "Reference coverage: " << ((double) reference.coverage)/reference.ref_genome.seq.size() << " " << ((double) reference.coverage_after_filter)/reference.ref_genome.seq.size() << std::endl; 
	}

	age_shared_count.resize(num_blocks);
	age_notshared_count.resize(num_blocks);
	age_shared_emp.resize(num_blocks);
	age_notshared_emp.resize(num_blocks);

	return num_blocks;

}


///////////////////////
void
maketmp_vcf(std::vector<std::string>& filename_chr, std::vector<std::string>& filename_mut, std::vector<std::string>& filename_target, std::vector<std::string>& filename_ref_genome, std::vector<std::string>& filename_mask, std::string& filename_output){

	int N_target, DAF, bp_target = 0, bp_mut = 0, num_used_snps = 0;
	bool target_flip = false;
	std::string ancestral, derived, ancestral_target, derived_target;
	Muts::iterator it_mut;

	int DAF_target, AAF_target;
	bool has_ref_genome = false; //new
	if(filename_ref_genome.size() > 0) has_ref_genome = true; //new
	bool has_mask = false;
	if(filename_mask.size() > 0) has_mask = true;

	FILE* fp = fopen(filename_output.c_str(), "wb");

	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

		Mutations mut;
		mut.Read(filename_mut[chr]);
		vcf_parser target(filename_target[chr]);
		fasta ref_genome; //new
		if(has_ref_genome) ref_genome.Read(filename_ref_genome[chr]); //new
		fasta mask;
		if(has_mask) mask.Read(filename_mask[chr]);

		bp_target = 0, bp_mut = 0;

		target.read_snp();
		bp_target = target.rec->pos + 1;
		target.extract_GT();
		N_target    = target.n * target.ploidy;

		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1){

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

					bool use = true;          
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;
					if(has_mask && bp_mut < mask.seq.size()){
						if(mask.seq[bp_mut-1] != 'P') use = false;
					}

					if(use){

						//check if mutation exists in ref
						if(bp_target < bp_mut){
							while(target.read_snp() == 0){
								bp_target = target.rec->pos + 1;
								target.extract_GT();  
								if(bp_target >= bp_mut) break;
							}
						}

						DAF_target = 0;
						if(bp_target == bp_mut){ //exists

							ancestral_target = target.rec->d.allele[0];
							derived_target   = target.rec->d.allele[1];

							bool accept = false, fixed_for_ref = false;
							target_flip = false;
							if(ancestral_target == derived && derived_target == ""){

								//possibility of sharing
								target.extract_GT();

								if(target.n*target.ploidy < 2){
									std::cerr << "Need at least one diploid genome." << std::endl;
									exit(1);
								}

								N_target = 0;

								bool biallelic = true;
								int k = 0;
								for(int i = 0; i < target.n; i++){
									for(int j = 0; j < target.ploidy; j++){
										DAF_target += bcf_gt_allele(target.gt[target.ploidy*i+j]);
										N_target++;
										k++;
										if(bcf_gt_allele(target.gt[target.ploidy*i+j]) > 1) biallelic = false;
									}
								}

								if(biallelic){
									if(DAF_target != 0){
										use = false;
									}else{
										DAF_target = N_target;
									}
								}else{
									use = false;
								}

							}else if( (ancestral_target == ancestral && derived_target == derived) || (ancestral_target == derived && derived_target == ancestral) ){

								target_flip = false;
								if( (ancestral_target == derived && derived_target == ancestral) ) target_flip = true;

								//possibility of sharing
								target.extract_GT();

								N_target = 0;

								bool biallelic = true;
								int k = 0;
								for(int i = 0; i < target.n; i++){
									for(int j = 0; j < target.ploidy; j++){
										DAF_target += bcf_gt_allele(target.gt[target.ploidy*i+j]);
										N_target++;
										k++;
										if(bcf_gt_allele(target.gt[target.ploidy*i+j]) > 1) biallelic = false;
									}
								}

								if(biallelic){
									if(target_flip){
										DAF_target = N_target - DAF_target;
									}
								}else{
									use = false;
								}
							}else{
								use = false;
							}
						}else{ //doesn't exist so use ref genome

							if(has_ref_genome){
								derived_target = ref_genome.seq[bp_mut-1];
								if(derived == derived_target){
									DAF_target = N_target;
								}else if(ancestral == derived_target){
									DAF_target = 0;
								}else{
									use = false;
								}
							}else{
								use = false;
							}

						}

						if(use){        
							//fprintf(fp, "%s %d %c %c %d %d", filename_chr[chr].c_str(), bp_mut, ancestral[0], derived[0], N_target-DAF_target, DAF_target);          
							int lchrom = strlen(filename_chr[chr].c_str());
							AAF_target = N_target - DAF_target;
							fwrite(&lchrom, sizeof(int), 1, fp);
							fwrite(filename_chr[chr].c_str(), sizeof(char), lchrom, fp);
							fwrite(&bp_mut, sizeof(int), 1, fp);
							fwrite(&ancestral[0], sizeof(char), 1, fp);
							fwrite(&derived[0], sizeof(char), 1, fp);
							fwrite(&AAF_target, sizeof(int), 1, fp);
							fwrite(&DAF_target, sizeof(int), 1, fp);
						}
					}

				}
			}
		}
	}
	fclose(fp);

}

void
maketmp_bam(std::string& params, std::vector<std::string>& filename_chr, std::vector<std::string>& filename_mut, std::string& filename_bam, std::vector<std::string>& filename_ref_genome, std::vector<std::string>& filename_mask, std::string& filename_output){

	int N_target, DAF, bp_target = 0, bp_mut = 0;
	bool target_flip = false;
	std::string ancestral, derived, ancestral_target, derived_target;
	Muts::iterator it_mut;

	int DAF_target, AAF_target;
	int bin_index;

	bam_parser target(filename_bam, params);
	FILE* fp = fopen(filename_output.c_str(), "wb");
	bool has_mask = false;
	if(filename_mask.size() > 0) has_mask = true;

	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr + 1 << " / " << filename_mut.size() << std::endl;

		target.assign_contig(filename_chr[chr], filename_ref_genome[chr]);

		fasta mask;
		if(has_mask) mask.Read(filename_mask[chr]);

		Mutations mut;
		mut.Read(filename_mut[chr]);

		bp_target = 0, bp_mut = 0;
		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1){

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

					bool use = true;
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(has_mask && bp_mut < mask.seq.size()){
						if(mask.seq[bp_mut-1] != 'P') use = false;
					}

					if(use){

						DAF_target = 0;
						AAF_target = 0;

						//check if mutation exists in target     
						int bp_target = bp_mut - 1;
						target.read_to_pos(bp_target);
						int num_alleles = 0;
						if(target.pos_of_entry[bp_target % target.num_entries] == bp_target){
							int num_reads = 0;
							for(int i = 0; i < 4; i++){
								num_reads += target.count_alleles[bp_target % target.num_entries][i];
								num_alleles += (target.count_alleles[bp_target % target.num_entries][i] > 0);
							}

							if(num_reads > 0){

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
									//if(!(num_alleles <= 2)) use = false;
									//if(AAF_target + DAF_target < 2) use = false;
								}else{
									use = false;
								}

							}else{
								use = false;
							}

						}else{
							use = false;
						}
					} 

					//
					//os << filename_chr[chr] << " " << bp_mut << " " << ancestral << " " << derived << " " << AAF_target << " " << DAF_target << "\n";
					//if(use) fprintf(fp, "%s %d %c %c %d %d", filename_chr[chr].c_str(), bp_mut, ancestral[0], derived[0], AAF_target, DAF_target);

					if(use){
						int lchrom = strlen(filename_chr[chr].c_str());
						fwrite(&lchrom, sizeof(int), 1, fp);
						fwrite(filename_chr[chr].c_str(), sizeof(char), lchrom, fp);
						fwrite(&bp_mut, sizeof(int), 1, fp);
						fwrite(&ancestral[0], sizeof(char), 1, fp);
						fwrite(&derived[0], sizeof(char), 1, fp);
						fwrite(&AAF_target, sizeof(int), 1, fp);
						fwrite(&DAF_target, sizeof(int), 1, fp);
					}

				}
			}
		}

		std::cerr << "Bam coverage: " << ((double) target.coverage)/target.ref_genome.seq.size() << " " << ((double) target.coverage_after_filter)/target.ref_genome.seq.size() << std::endl;
	}
	fclose(fp);

}

int
parse_tmptmp(std::vector<std::string>& filename_chr, std::vector<std::string>& filename_mut, std::string& filename_target, std::string& filename_ref, std::vector<std::string>& filename_target_mask, std::vector<std::string>& filename_reference_mask, double age, double ref_age, double C, std::mt19937& rng, int num_bases_per_block, std::vector<std::vector<double>>& age_shared_count, std::vector<std::vector<double>>& age_notshared_count, std::vector<std::vector<double>>& age_shared_emp, std::vector<std::vector<double>>& age_notshared_emp){

	std::uniform_real_distribution<double> dist_unif(0,1);
	int num_age_bins = ((int) (log(1e8) * C))+1;
	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0, num_used_snps2 = 0;
	bool ref_flip = false, target_flip = false;
	std::string ancestral, derived;
	char ancestral_ref, derived_ref, ancestral_target, derived_target;
	Muts::iterator it_mut;

	float num_samples = 100;
	char chrom_target[1024], chrom_ref[1024];
	int DAF_target, AAF_target, DAF_ref, AAF_ref, lchrom;
	int bin_index;

	FILE* fp_target = fopen(filename_target.c_str(), "rb");
	FILE* fp_ref    = fopen(filename_ref.c_str(), "rb");

	bool has_tar_mask = false, has_ref_mask = false;
	if(filename_target_mask.size() > 0) has_tar_mask = true;
	if(filename_reference_mask.size() > 0) has_ref_mask = true;

	std::vector<std::vector<double>>::iterator it_age_shared_count    = age_shared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();
	std::vector<std::vector<double>>::iterator it_age_shared_emp      = age_shared_emp.begin();
	std::vector<std::vector<double>>::iterator it_age_notshared_emp   = age_notshared_emp.begin();
	int current_block_base = 0;
	int num_blocks = 0;

	for(int chr = 0; chr < filename_mut.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << filename_mut.size() << std::endl;

		Mutations mut;
		mut.Read(filename_mut[chr]);

		fasta tar_mask, ref_mask;
		if(has_tar_mask) tar_mask.Read(filename_target_mask[chr]);
		if(has_ref_mask) ref_mask.Read(filename_reference_mask[chr]);

		current_block_base = 0;
		bp_mut = 0;

		while( strcmp(chrom_ref, filename_chr[chr].c_str()) != 0 ){
			if(fread(&lchrom, sizeof(int), 1, fp_ref) != 1) break;
			fread(chrom_ref, sizeof(char), lchrom, fp_ref);
			chrom_ref[lchrom] = '\0';
			fread(&bp_ref, sizeof(int), 1, fp_ref);
			fread(&ancestral_ref, sizeof(char), 1, fp_ref);
			fread(&derived_ref, sizeof(char), 1, fp_ref);
			fread(&AAF_ref, sizeof(int), 1, fp_ref);
			fread(&DAF_ref, sizeof(int), 1, fp_ref);
		}

		while( strcmp(chrom_target, filename_chr[chr].c_str()) != 0 ){
			if(fread(&lchrom, sizeof(int), 1, fp_target) != 1) break;
			fread(chrom_target, sizeof(char), lchrom, fp_target);
			chrom_target[lchrom] = '\0';
			fread(&bp_target, sizeof(int), 1, fp_target);
			fread(&ancestral_target, sizeof(char), 1, fp_target);
			fread(&derived_target, sizeof(char), 1, fp_target);
			fread(&AAF_target, sizeof(int), 1, fp_target);
			fread(&DAF_target, sizeof(int), 1, fp_target);
		}



		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0 && (*it_mut).age_end >= age){

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

					bool use = true;
					if(has_tar_mask && bp_mut < tar_mask.seq.size()){
						if(tar_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(has_ref_mask && bp_mut < ref_mask.seq.size()){
						if(ref_mask.seq[bp_mut-1] != 'P') use = false;
					}
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;

					//std::cerr << filename_chr[chr] << " " << chrom_ref << " " << chrom_target << std::endl;
					//std::cerr << bp_mut << " " << bp_target << " " << bp_ref << std::endl;

					if(use){
						DAF_ref    = 0;
						AAF_ref    = 0;
						while( (strcmp(chrom_ref, filename_chr[chr].c_str()) == 0 && bp_ref < bp_mut) ){
							if(fread(&lchrom, sizeof(int), 1, fp_ref) != 1) break;
							fread(chrom_ref, sizeof(char), lchrom, fp_ref);
							chrom_ref[lchrom] = '\0';
							fread(&bp_ref, sizeof(int), 1, fp_ref);
							fread(&ancestral_ref, sizeof(char), 1, fp_ref);
							fread(&derived_ref, sizeof(char), 1, fp_ref);
							fread(&AAF_ref, sizeof(int), 1, fp_ref);
							fread(&DAF_ref, sizeof(int), 1, fp_ref);
						}
						if( strcmp(chrom_ref, filename_chr[chr].c_str()) != 0 || bp_ref != bp_mut || ancestral_ref != ancestral[0] || derived_ref != derived[0]){
							use = false;
						}
					} 
					if(DAF_ref == 0) use = false;
					N_ref = DAF_ref + AAF_ref;

					if(use){
						DAF_target = 0;
						AAF_target = 0;
						while( (strcmp(chrom_target, filename_chr[chr].c_str()) == 0 && bp_target < bp_mut) ){
							if(fread(&lchrom, sizeof(int), 1, fp_target) != 1) break;
							fread(chrom_target, sizeof(char), lchrom, fp_target);
							chrom_target[lchrom] = '\0';
							fread(&bp_target, sizeof(int), 1, fp_target);
							fread(&ancestral_target, sizeof(char), 1, fp_target);
							fread(&derived_target, sizeof(char), 1, fp_target);
							fread(&AAF_target, sizeof(int), 1, fp_target);
							fread(&DAF_target, sizeof(int), 1, fp_target);
						}
						if( strcmp(chrom_target, filename_chr[chr].c_str()) != 0 || bp_target != bp_mut || ancestral_target != ancestral[0] || derived_target != derived[0]){
							use = false;
						}
					} 
					N_target = DAF_target + AAF_target;

					//populate age_shared_count and age_notshared_count
					if(use){

						double age_begin = (*it_mut).age_begin;
						if(age_begin < ref_age) age_begin = ref_age;

						while(current_block_base + num_bases_per_block < bp_mut){
							current_block_base += num_bases_per_block;
							it_age_shared_emp++;
							it_age_notshared_emp++;
							it_age_shared_count++;
							it_age_notshared_count++;
							num_blocks++;
						}

						bool skip = false;
						if(age_begin <= age){

							if(1){
								double age_begin2 = age_begin;
								age_begin2 = 0.0;
								int bin_index1 = std::max(0, (int)std::round(log(10*age_begin2)*C)+1);
								int bin_index2 = std::max(0, (int)std::round(log(10*(*it_mut).age_end)*C)+1);
								bin_index  = bin_index1 * num_age_bins + bin_index2;
								(*it_age_shared_emp)[bin_index]    += DAF_target * DAF_ref/((double)N_ref);
								(*it_age_notshared_emp)[bin_index] += AAF_target * DAF_ref/((double)N_ref);

								if(1){
									int j = 0;
									while(j < num_samples){
										skip = false;
										double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
										//if(sampled_age < age) skip = true;
										if(sampled_age < age) sampled_age = age;
										int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
										bin_index = bin_index_age;// * num_age_bins + bin_index_age;

										if(!skip){
											(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
											//j++;
										}
										j++;
									}
								}
							}

						}else{

							int j = 0;
							while(j < num_samples){
								skip = false;
								double sampled_age = dist_unif(rng) * ((*it_mut).age_end - age_begin) + age_begin;
								if(sampled_age < age) skip = true;
								int bin_index_age = std::max(0, (int)std::round(log(10*sampled_age)*C)+1);
								bin_index = bin_index_age;// * num_age_bins + bin_index_age;

								//if(bin_index >= (*it_age_shared_count).size()) std::cerr << (*it_mut).age_begin << " " << (*it_mut).age_end << " " << sampled_age << std::endl;
								//assert(bin_index < (*it_age_shared_count).size());
								if(bin_index >= (*it_age_shared_count).size()) skip = true;
								if(!skip){
									(*it_age_shared_count)[bin_index]    += DAF_target * DAF_ref/((double)N_ref * num_samples);
									(*it_age_notshared_count)[bin_index] += AAF_target * DAF_ref/((double)N_ref * num_samples);
									j++;
								}
							}

						}


					}

				}
			}
		}

		it_age_shared_count++;
		it_age_notshared_count++;
		it_age_shared_emp++;
		it_age_notshared_emp++;
		num_blocks++;
	}
	fclose(fp_target);
	fclose(fp_ref);
	age_shared_count.resize(num_blocks);
	age_notshared_count.resize(num_blocks);
	age_shared_emp.resize(num_blocks);
	age_notshared_emp.resize(num_blocks);

	return num_blocks;

}

//////////////////////
void
calc_depth(cxxopts::Options& options){

	//Program options

	bool help = false;
	if(!options.count("target_bam") || !options.count("ref_genome") || !options.count("chr") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: target_bam, ref_genome, chr, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Calculate target bam." << std::endl;
		exit(0);
	}  

	std::string params = "20,30,10";
	if(options.count("filters")){
		params = options["filters"].as<std::string>();
	}

	igzstream is_chr(options["chr"].as<std::string>());
	if(is_chr.fail()){
		std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
	}
	std::string line;
	std::vector<std::string> name_chr, filename_ref_genome;
	while(getline(is_chr, line)){
		name_chr.push_back(line);
		filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
	}
	is_chr.close();

	bam_parser target(options["target_bam"].as<std::string>(), params);

	double cov = 0.0, cov_filtered = 0.0, genome_length = 0.0;

	for(int chr = 0; chr < name_chr.size(); chr++){
		target.assign_contig(name_chr[chr], filename_ref_genome[chr]);
		target.read_to_pos(target.ref_genome.seq.size());

		cov           += target.coverage;
		cov_filtered  += target.coverage_after_filter;
		genome_length += target.ref_genome.seq.size();
	}

	std::ofstream os(options["output"].as<std::string>());
	os << cov/genome_length << " " << cov_filtered/genome_length << std::endl;
	os.close();

}

void
make_tmp(cxxopts::Options& options){

	//Program options

	bool help = false;
	if(!options.count("mut") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: mut, ref_genome, output, either of target_bcf or target_bam. Optional: target_mask, filters." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Calculate coalescence rates for sample." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating Colate tmp input file for " << options["target_tmp"].as<std::string>() << ".." << std::endl;

	std::string params = "20,30,10";
	if(options.count("filters")){
		params = options["filters"].as<std::string>();
	}

	std::string line;
	////////////////////////////////////////

	std::vector<std::string> filename_mut, filename_target, filename_ref_genome, filename_mask, name_chr;

	if(options.count("target_bcf") > 0){
		if(options.count("chr") > 0){

			igzstream is_chr(options["chr"].as<std::string>());
			if(is_chr.fail()){
				std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
			}

			if(options.count("target_mask") == 0){
        std::cerr << "Warning: no mask file given. Will assume that all sites not in BCF are homozygous REF." << std::endl;
			}

			while(getline(is_chr, line)){
				name_chr.push_back(line);
				filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
				filename_target.push_back(options["target_bcf"].as<std::string>() + "_chr" + line + ".bcf");
				filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
				if(options.count("target_mask")) filename_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
			}
			is_chr.close();

		}else{

			name_chr.push_back("");
			filename_mut.push_back(options["mut"].as<std::string>());
			filename_target.push_back(options["target_bcf"].as<std::string>());
			filename_ref_genome.push_back(options["ref_genome"].as<std::string>());
			if(options.count("target_mask")) filename_mask.push_back(options["target_mask"].as<std::string>());
		}

		std::string filename_output = options["output"].as<std::string>() + ".colate.in";
		maketmp_vcf(name_chr, filename_mut, filename_target, filename_ref_genome, filename_mask, filename_output);
	}

	if(options.count("target_bam") > 0){
		if(options.count("chr") > 0){

			igzstream is_chr(options["chr"].as<std::string>());
			if(is_chr.fail()){
				std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
			}
			while(getline(is_chr, line)){
				name_chr.push_back(line);
				filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
				filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
				if(options.count("target_mask")) filename_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
			}
			is_chr.close();

		}else{

			name_chr.push_back("");
			filename_mut.push_back(options["mut"].as<std::string>());
			filename_ref_genome.push_back(options["ref_genome"].as<std::string>());
			if(options.count("target_mask")) filename_mask.push_back(options["target_mask"].as<std::string>());

		}

		std::string filename_bam    = options["target_bam"].as<std::string>();
		std::string filename_output = options["output"].as<std::string>() + ".colate.in";
		maketmp_bam(params, name_chr, filename_mut, filename_bam, filename_ref_genome, filename_mask, filename_output);
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
mut(cxxopts::Options& options){

	//Program options

	bool help = false;
	if(!options.count("mut") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: mut, bins, output. Optional: target_tmp, reference_tmp, target_bcf, reference_bcf, target_bam, reference_bam, ref_genome, target_age, reference_age, target_mask, reference_mask, coal, num_bootstrap, filters." << std::endl;
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
	//if(options.count("regularise") > 0){ 
	//  regularise = options["regularise"].as<int>();
	//}

	std::string params = "20,30,10";
	if(options.count("filters")){
		params = options["filters"].as<std::string>();
	}

	bool is_ancient = false;
	double age;
	double log_10 = std::log(10);

	double target_age = 0;
	double ref_age    = 0;
	if(options.count("target_age") > 0) target_age = std::stof(options["target_age"].as<std::string>());
	if(options.count("reference_age") > 0) ref_age    = std::stof(options["reference_age"].as<std::string>());
	assert(target_age >= 0.0);
	assert(ref_age >= 0.0);

	double years_per_gen = 28.0;
	if(options.count("years_per_gen") > 0){
		years_per_gen = options["years_per_gen"].as<float>();
	}

	age = std::max(target_age, ref_age)/years_per_gen;
	ref_age /= years_per_gen;
	std::cerr << age << std::endl; 
	if(age > 0.0) is_ancient = true;

	std::string line;
	////////////////////////////////////////

	//formula for converting age to int: log(age)*C, Ne = haploid population size
	double C = 10;
	int num_age_bins = ((int) (log(1e8) * C))+1;
	std::cerr << "num_bins: " << num_age_bins << std::endl;
	std::vector<double> age_bin(num_age_bins, 0.0);
	std::vector<double>::iterator it_age_bin = age_bin.begin();
	int bin = 0;
	*it_age_bin = 0.0;
	it_age_bin++;
	for(; bin < num_age_bins-1; bin++){
		*it_age_bin =  exp(bin/C)/10.0;
		it_age_bin++;
	}

	int num_bases_per_block = 10e6;
	int num_blocks = 500;

	std::vector<std::vector<double>> age_shared_count_block(num_blocks), age_notshared_count_block(num_blocks);
	std::vector<std::vector<double>> age_shared_emp_block(num_blocks), age_notshared_emp_block(num_blocks);
	for(int i = 0; i < num_blocks; i++){
		age_shared_count_block[i].resize(num_age_bins);
		age_notshared_count_block[i].resize(num_age_bins);
		age_shared_emp_block[i].resize(num_age_bins*num_age_bins);
		age_notshared_emp_block[i].resize(num_age_bins*num_age_bins);
		std::fill(age_shared_count_block[i].begin(), age_shared_count_block[i].end(), 0.0);
		std::fill(age_notshared_count_block[i].begin(), age_notshared_count_block[i].end(), 0.0);
		std::fill(age_shared_emp_block[i].begin(), age_shared_emp_block[i].end(), 0.0);
		std::fill(age_notshared_emp_block[i].begin(), age_notshared_emp_block[i].end(), 0.0);
	}

	///////////

	std::mt19937 rng;
	int seed = std::time(0) + getpid();
	if(options.count("seed") > 0){
		seed = options["seed"].as<int>();
	}
	rng.seed(seed);

	int num_bootstrap = 1;
	if(options.count("num_bootstraps") > 0) num_bootstrap = options["num_bootstraps"].as<int>();
	std::vector<std::vector<double>> age_shared_count(num_bootstrap), age_notshared_count(num_bootstrap);
	std::vector<double> age_shared_emp(num_age_bins*num_age_bins), age_notshared_emp(num_age_bins*num_age_bins);

	igzstream is_mat(options["output"].as<std::string>() + ".colate_mat");
	if(is_mat.fail()){

		///////////
		std::vector<std::string> filename_mut, filename_target, filename_ref, filename_ref_genome, filename_target_mask, filename_reference_mask, name_chr;

		if(options.count("target_bcf") && options.count("reference_bcf")){
			if(options.count("chr") > 0){

				igzstream is_chr(options["chr"].as<std::string>());
				if(is_chr.fail()){
					std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
				}
				while(getline(is_chr, line)){
					filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
					filename_target.push_back(options["target_bcf"].as<std::string>() + "_chr" + line + ".bcf");
					filename_ref.push_back(options["reference_bcf"].as<std::string>() + "_chr" + line + ".bcf");
					if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("ref_genome") > 0) filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
				}
				is_chr.close();

			}else{
				filename_mut.push_back(options["mut"].as<std::string>());
				filename_target.push_back(options["target_bcf"].as<std::string>());
				filename_ref.push_back(options["reference_bcf"].as<std::string>());
				if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>());
				if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>());
				if(options.count("ref_genome") > 0) filename_ref_genome.push_back(options["ref_genome"].as<std::string>());
			}
			num_blocks = parse_vcfvcf(filename_mut, filename_target, filename_ref, filename_target_mask, filename_reference_mask, filename_ref_genome, age, ref_age, C, rng, num_bases_per_block, age_shared_count_block, age_notshared_count_block, age_shared_emp_block, age_notshared_emp_block);

		}else if(options.count("target_bcf")){
			if(options.count("chr") > 0){

				igzstream is_chr(options["chr"].as<std::string>());
				if(is_chr.fail()){
					std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
				}
				while(getline(is_chr, line)){
					filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
					filename_target.push_back(options["target_bcf"].as<std::string>() + "_chr" + line + ".bcf");
					if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("ref_genome") > 0) filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
				}
				is_chr.close();

			}else{
				filename_mut.push_back(options["mut"].as<std::string>());
				filename_target.push_back(options["target_bcf"].as<std::string>());
				if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>());
				if(options.count("ref_genome") > 0) filename_ref_genome.push_back(options["ref_genome"].as<std::string>());
			}
			num_blocks = parse_vcf(filename_mut, filename_target, filename_target_mask, filename_ref_genome, age, ref_age, C, rng, num_bases_per_block, age_shared_count_block, age_notshared_count_block, age_shared_emp_block, age_notshared_emp_block);

		}else if(options.count("target_bam") && options.count("reference_bcf")){

			filename_target.push_back(options["target_bam"].as<std::string>() + ".bam");
			//check order of contigs
			//check coverage for each contig


			if(options.count("chr") > 0){

				igzstream is_chr(options["chr"].as<std::string>());
				if(is_chr.fail()){
					std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
				}
				while(getline(is_chr, line)){
					name_chr.push_back(line);
					filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
					filename_ref.push_back(options["reference_bcf"].as<std::string>() + "_chr" + line + ".bcf");
					filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>() + "_chr" + line + ".fa");
				}
				is_chr.close();

			}else{
				name_chr.push_back("");
				filename_mut.push_back(options["mut"].as<std::string>());
				filename_ref.push_back(options["reference_bcf"].as<std::string>());
				filename_ref_genome.push_back(options["ref_genome"].as<std::string>());
				if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>());
				if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>());
			}
			num_blocks = parse_onebamvcf(params, name_chr, filename_mut, filename_target, filename_ref, filename_target_mask, filename_reference_mask, filename_ref_genome, age, ref_age, C, rng, num_bases_per_block, age_shared_count_block, age_notshared_count_block, age_shared_emp_block, age_notshared_emp_block);

		}else if(options.count("target_bam") && options.count("reference_bam")){

			filename_target.push_back(options["target_bam"].as<std::string>());
			filename_ref.push_back(options["reference_bam"].as<std::string>());
			//check order of contigs
			//check coverage for each contig

			if(options.count("chr") > 0){

				igzstream is_chr(options["chr"].as<std::string>());
				if(is_chr.fail()){
					std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
				}
				while(getline(is_chr, line)){
					name_chr.push_back(line);
					filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
					//filename_target.push_back(options["target_bam"].as<std::string>() + "_chr" + line + ".bam");
					filename_ref_genome.push_back(options["ref_genome"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>() + "_chr" + line + ".fa");
				}
				is_chr.close();

			}else{
				name_chr.push_back("");
				filename_mut.push_back(options["mut"].as<std::string>());
				filename_ref_genome.push_back(options["ref_genome"].as<std::string>());
				if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>());
				if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>());
			}
			num_blocks = parse_onebambam(params, name_chr, filename_mut, filename_target, filename_ref, filename_target_mask, filename_reference_mask, filename_ref_genome, age, ref_age, C, rng, num_bases_per_block, age_shared_count_block, age_notshared_count_block, age_shared_emp_block, age_notshared_emp_block);

		}else if(options.count("target_tmp") && options.count("reference_tmp")){

			filename_target.push_back(options["target_tmp"].as<std::string>());
			filename_ref.push_back(options["reference_tmp"].as<std::string>());
			//check order of contigs
			//check coverage for each contig

			if(options.count("chr") > 0){

				igzstream is_chr(options["chr"].as<std::string>());
				if(is_chr.fail()){
					std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
				}
				while(getline(is_chr, line)){
					name_chr.push_back(line);
					filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
					if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>() + "_chr" + line + ".fa");
					if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>() + "_chr" + line + ".fa");
				}
				is_chr.close();

			}else{
				name_chr.push_back("");
				filename_mut.push_back(options["mut"].as<std::string>());
				if(options.count("target_mask") > 0) filename_target_mask.push_back(options["target_mask"].as<std::string>());
				if(options.count("reference_mask") > 0) filename_reference_mask.push_back(options["reference_mask"].as<std::string>());
			}
			num_blocks = parse_tmptmp(name_chr, filename_mut, filename_target[0], filename_ref[0], filename_target_mask, filename_reference_mask, age, ref_age, C, rng, num_bases_per_block, age_shared_count_block, age_notshared_count_block, age_shared_emp_block, age_notshared_emp_block);

		}

		std::cerr << "Number of blocks: " << num_blocks << std::endl;

		bool tmp_file = true;
		if(options.count("target_tmp") > 0 && options.count("reference_tmp") > 0) tmp_file = false;

		for(int i = 0; i < num_bootstrap; i++){
			age_shared_count[i].resize(num_age_bins);
			age_notshared_count[i].resize(num_age_bins);
		}
		std::uniform_int_distribution<int> dist_blocks(0,num_blocks-1);
		std::vector<double> blocks(num_blocks);
		std::vector<std::vector<double>>::iterator it_age_shared_count = age_shared_count.begin();
		std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();
		std::vector<double>::iterator it1, it2;
		//print matrix
		std::ofstream os_mat;
		if(tmp_file){
			os_mat.open(options["output"].as<std::string>() + ".colate_mat");  
			for(int bin = 0; bin < num_age_bins; bin++){
				os_mat << age_bin[bin] << " "; 
			}
			os_mat << "\n";
		}
		for(int i = 0; i < num_bootstrap; i++){
			std::fill(age_shared_count[i].begin(), age_shared_count[i].end(), 0.0);
			std::fill(age_notshared_count[i].begin(), age_notshared_count[i].end(), 0.0);
			std::fill(age_shared_emp.begin(), age_shared_emp.end(), 0.0);
			std::fill(age_notshared_emp.begin(), age_notshared_emp.end(), 0.0);

			if(num_bootstrap == 1){
				std::fill(blocks.begin(), blocks.end(), 1.0);
			}else{
				std::fill(blocks.begin(), blocks.end(), 0.0);
				for(int j = 0; j < num_blocks; j++){
					blocks[dist_blocks(rng)] += 1.0;
				}
			}
			for(int j = 0; j < num_blocks; j++){
				if(blocks[j] > 0.0){
					it1 = (*it_age_shared_count).begin();
					it2 = (age_shared_count_block[j]).begin();
					for(; it1 != (*it_age_shared_count).end();){
						(*it1) += blocks[j]*(*it2);
						it1++;
						it2++;
					}
					it1 = (*it_age_notshared_count).begin();
					it2 = (age_notshared_count_block[j]).begin();
					for(; it1 != (*it_age_notshared_count).end();){
						(*it1) += blocks[j]*(*it2);
						it1++;
						it2++;
					}

					it1 = (age_shared_emp).begin();
					it2 = (age_shared_emp_block[j]).begin();
					for(; it1 != (age_shared_emp).end();){
						(*it1) += blocks[j]*(*it2);
						it1++;
						it2++;
					}
					it1 = (age_notshared_emp).begin();
					it2 = (age_notshared_emp_block[j]).begin();
					for(; it1 != (age_notshared_emp).end();){
						(*it1) += blocks[j]*(*it2);
						it1++;
						it2++;
					}
				}
			}

			std::vector<double> F(num_age_bins, 0.0), G(num_age_bins, 0.0);
			double fcount, gcount;
			int bin = 0, bin_start = 0;
			while(age_bin[bin] <= age) bin++;
			bin_start = bin;
			for(int bin1 = 0; bin1 < 1; bin1++){

				double lower_age = age_bin[bin1];
				if(bin1 == 0) lower_age = age_bin[bin_start-1];
				std::fill(F.begin(), F.end(), 0.0);
				std::fill(G.begin(), G.end(), 0.0);
				fcount = 0.0;
				gcount = 0.0;

				bin = bin_start;     
				for(; bin < num_age_bins; bin++){ 
					fcount += age_shared_emp[bin1*num_age_bins+bin];
					gcount += age_notshared_emp[bin1*num_age_bins+bin];
					if(age_shared_emp[bin1*num_age_bins+bin] > 0){
						F[bin] = age_shared_emp[bin1*num_age_bins+bin]/(age_shared_emp[bin1*num_age_bins+bin] + age_notshared_emp[bin1*num_age_bins+bin]);
					}
					if(age_notshared_emp[bin1*num_age_bins+bin] > 0){
						G[bin] = age_notshared_emp[bin1*num_age_bins+bin];//(age_shared_emp[bin1*num_age_bins+bin] + age_notshared_emp[bin1*num_age_bins+bin]);
					}
					//std::cerr << (age_shared_emp[bin1*num_age_bins+bin] + age_notshared_emp[bin1*num_age_bins+bin])/(age_bin[bin+1] - age_bin[bin]) << " ";
				}

				if(1){
					bin = bin_start;
					for(; bin < num_age_bins; bin++){
						F[bin-1] *= (age_bin[bin] - lower_age);
						//G[bin-1] *= (age_bin[bin] - lower_age);
						lower_age = age_bin[bin];
					}
				}

				double normf = 0.0, normg = 0.0;
				for(bin = 0; bin < num_age_bins; bin++){
					normf += F[bin];
					normg += G[bin];
				}

				if(1){
					for(bin = 0; bin < num_age_bins; bin++){
						//F[bin]  /= F[num_age_bins-1];
						F[bin] /= normf;
						F[bin] *= fcount;
						int bin_index = bin;//*num_age_bins+bin;
						(*it_age_shared_count)[bin_index] += std::max(0.0, F[bin]);
					}
				}
				if(0){
					for(bin = 0; bin < num_age_bins; bin++){
						G[bin] /= normg;
						G[bin] *= gcount;
						int bin_index = bin;//*num_age_bins+bin;
						(*it_age_notshared_count)[bin_index] += std::max(0.0, G[bin]);
					}
				}
			}

			double norm = 1e3;
			if(tmp_file){
				for(int bin = 0; bin < num_age_bins; bin++){
					(*it_age_shared_count)[bin] /= norm;
					os_mat << (*it_age_shared_count)[bin] << " ";
				}
				os_mat << "\n";
				for(int bin = 0; bin < num_age_bins; bin++){
					(*it_age_notshared_count)[bin] /= norm;
					os_mat << (*it_age_notshared_count)[bin] << " ";
				}
				os_mat << "\n";
			}
			it_age_shared_count++;
			it_age_notshared_count++;

		}
		if(tmp_file) os_mat.close();
	}else{
		std::cerr << "Loading precomputed file " << options["output"].as<std::string>() << ".colate_mat" << std::endl;

		for(int i = 0; i < num_bootstrap; i++){
			age_shared_count[i].resize(num_age_bins);
			age_notshared_count[i].resize(num_age_bins);
		}
		std::vector<std::vector<double>>::iterator it_age_shared_count = age_shared_count.begin();
		std::vector<std::vector<double>>::iterator it_age_notshared_count = age_notshared_count.begin();

		for(int bin = 0; bin < num_age_bins; bin++){
			is_mat >> age_bin[bin]; 
		}
		for(; it_age_shared_count != age_shared_count.end();){

			std::fill((*it_age_shared_count).begin(), (*it_age_shared_count).end(), 0.0);
			std::fill((*it_age_notshared_count).begin(), (*it_age_notshared_count).end(), 0.0);

			for(int bin = 0; bin < num_age_bins; bin++){
				is_mat >> (*it_age_shared_count)[bin];
			}
			for(int bin = 0; bin < num_age_bins; bin++){
				is_mat >> (*it_age_notshared_count)[bin];
			}
			it_age_shared_count++;
			it_age_notshared_count++;
		}
		is_mat.close();
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
		std::string tmp = "";
		int ep = 0;
		for(int i = 0; i < line.size(); i++){
			if(line[i] == ' ' || line[i] == '\t'){
				if(ep == 1 && age < stof(tmp) && age != 0.0){
					epochs.push_back(age);
					ep++;
				}
				if(ep != 1 || age == 0.0){
					epochs.push_back(stof(tmp));
					ep++;
				}
				tmp = "";
			}else{
				tmp += line[i];
			}
		}

		if(tmp != ""){
			if(ep == 1 && age < stof(tmp) && age != 0.0){
				epochs.push_back(age);
				ep++;
			}
			if(ep != 1 || age == 0.0){
				epochs.push_back(stof(tmp));
				ep++;
			}
		}
		if(age != 0.0) assert(ep > 2);
		if(age == 0.0) assert(ep > 1);

		num_epochs = epochs.size();

		assert(epochs[0] == 0);
		for(int e = 1; e < num_epochs; e++){
			assert(epochs[e] > epochs[e-1]);
		}

	}else{ 

		double count = 0.0;
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
					if(epoch_boundary - log_age < 0.1*epoch_step) epoch_boundary += epoch_step;
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
	std::vector<double> coal_rates(num_epochs, initial_coal_rate), coal_rates_init(num_epochs, initial_coal_rate), coal_rates_num(num_epochs, 0.0), coal_rates_denom(num_epochs, 0.0);
	if(options.count("coal") > 0){
		double dummy;
		is >> dummy >> dummy;
		for(int e = 0; e < num_epochs; e++){
			is >> coal_rates_init[e];
			std::cerr << coal_rates_init[e] << " "; 
		}
		std::cerr << std::endl;
	}


	////////////////////////////////////////

	std::cerr << "Maximising likelihood using EM.. " << std::endl;


	//coal_EM EM(epochs, coal_rates);

	int max_iter = 100000;
	//max_iter = 10000;

	//std::ofstream os_log(options["output"].as<std::string>() + ".log");
	std::ofstream os(options["output"].as<std::string>() + ".coal");
	os << "0\n";
	for(int e = 0; e < num_epochs; e++){
		os << epochs[e] << " ";
	}
	os << std::endl;

	for(int i = 0; i < num_bootstrap; i++){

		double gamma = 1;
		std::vector<double> gradient(num_epochs, 0.0), gradient_prev(num_epochs, 0.0), coal_rates_prev(num_epochs, 0.0), coal_rates_nesterov(num_epochs, 0.0);
		coal_rates = coal_rates_init;
		std::fill(coal_rates_num.begin(), coal_rates_num.end(), 0.0);
		std::fill(coal_rates_denom.begin(), coal_rates_denom.end(), 0.0);
		coal_rates_nesterov = coal_rates;

		int perc = -1;
		double log_likelihood = log(0.0), prev_log_likelihood = log(0.0);
		for(int iter = 0; iter < max_iter; iter++){

			if( (int) (((double)iter)/max_iter * 100.0) > perc ){
				perc = (int) (((double)iter)/max_iter * 100.0);
				//std::cerr << "[" << perc << "%]\r";
				std::cerr << "Iteration: " << iter << "\r";
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
			for(int bin = 0; bin < num_age_bins; bin++){

				if(age_shared_count[i][bin] > 0){ 
					count = age_shared_count[i][bin];
					double logl = EM.EM_shared(age_bin[bin], age_bin[bin], num, denom);
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
				if(age_notshared_count[i][bin] > 0){ 
					count = age_notshared_count[i][bin];
					double logl = EM.EM_notshared(age_bin[bin], age_bin[bin], num, denom);
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
						if(coal_rates[e] < 5e-9){
							coal_rates[e] = 5e-9;
						}
						if(coal_rates_nesterov[e] < 5e-9){
							coal_rates_nesterov[e] = 5e-9;
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

			//os_log << log_likelihood << " " << log_likelihood - prev_log_likelihood << std::endl;
			std::fill(coal_rates_num.begin(), coal_rates_num.end(), 0.0);
			std::fill(coal_rates_denom.begin(), coal_rates_denom.end(), 0.0);

			if(log_likelihood/prev_log_likelihood > 1.0 - 1e-10 && iter > 1000){
				std::cerr << "Iteration: " << iter << std::endl;
				break;
			}

		}
		//os_log.close();

		os << "0 " << i << " ";
		for(int e = 0; e < num_epochs; e++){
			os << coal_rates[e] << " ";
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
preprocess_mut(cxxopts::Options& options){

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("reference_bcf") || !options.count("ref_genome") || !options.count("anc_genome") || !options.count("mask") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, reference_bcf, ref_genome, anc_genome, mask, output. Optional: outgroup_tmrca, years_per_gen" << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Make mut file including fixed mutations." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Adding fixed mutations to mut file.." << std::endl;

	float years_per_gen = 28.0;
	if(options.count("years_per_gen")){
		years_per_gen = options["years_per_gen"].as<float>();
	}

	double outgroup_age = 10e6/years_per_gen;
	if(options.count("outgroup_tmrca")){
		outgroup_age = options["outgroup_tmrca"].as<float>()/years_per_gen;
	}

	fasta ref_genome, anc_genome, mask;
	ref_genome.Read(options["ref_genome"].as<std::string>());
	anc_genome.Read(options["anc_genome"].as<std::string>());
	mask.Read(options["mask"].as<std::string>());

	vcf_parser vcf(options["reference_bcf"].as<std::string>());

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
	int tmp = (*it_mut).pos;
	while(vcf.read_snp() == 0){

		if(mut_combined.info.size() - snp_comb < L_ref) mut_combined.info.resize(mut_combined.info.size() + L_ref);

		vcf.extract_GT();
		int N = vcf.ploidy*vcf.n;
		if(0){
			if(N != ancmut.NumTips()){
				std::cerr << N << " " << ancmut.NumTips() << std::endl;
				std::cerr << "Number of samples in vcf and anc/mut don't match" << std::endl;
				exit(1);
			}
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

		if(1){
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
		} 

		if(biallelic){

			//check if snp exists in mut
			if(snp_ref < L_ref){
				while((*it_mut).pos < bp+1){
					//tmp = (*it_mut).pos;
					num_bases_snp_persists = ancmut.NextSNP(mtr, it_mut);
					snp_ref++;
					if(snp_ref == L_ref) break;
				}
				//if((*it_mut).pos == tmp) biallelic = false;
				//tmp = (*it_mut).pos;
			}

			if(tree_count < (*it_mut).tree){
				tree_count = (*it_mut).tree;
				mtr.tree.GetCoordinates(coords);
				tmrca = coords[root];
			} 

			//if(bp == 17071275) std::cerr << (*it_mut).pos << " " << (*it_mut).branch.size() << " " << (*it_mut).flipped << " " << DAF << " " << biallelic << std::endl;

			if((*it_mut).pos == bp+1 && DAF > 0 && DAF < N && biallelic){

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

	std::cerr << "Outputting " << snp_comb << " SNPs." << std::endl;
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
print_tmp(cxxopts::Options& options){

	//Program options

	bool help = false;
	if(!options.count("target_tmp") || !options.count("mut") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: target_tmp, mut, output. Optional: chr." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Print tmp." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating Colate tmp input file for " << options["target_tmp"].as<std::string>() << ".." << std::endl;

	std::string params = "20,30,10";
	if(options.count("filters")){
		params = options["filters"].as<std::string>();
	}

	std::string line;
	////////////////////////////////////////

	std::vector<std::string> name_chr, filename_mut;

	if(options.count("chr") > 0){

		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		while(getline(is_chr, line)){
			name_chr.push_back(line);
			filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
		}
		is_chr.close();

	}else{

		name_chr.push_back("");
		filename_mut.push_back(options["mut"].as<std::string>());

	}

	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0, num_used_snps2 = 0;
	std::string ancestral, derived;
	char ancestral_ref, derived_ref, ancestral_target, derived_target;
	Muts::iterator it_mut;

	char chrom_target[1024], chrom_ref[1024];
	int DAF_target, AAF_target, DAF_ref, AAF_ref, lchrom;

	FILE* fp_target = fopen(options["target_tmp"].as<std::string>().c_str(), "rb");

  std::ofstream os(options["output"].as<std::string>());

	for(int chr = 0; chr < name_chr.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << name_chr.size() << std::endl;

		Mutations mut;
		mut.Read(filename_mut[chr]);

		while( strcmp(chrom_target, name_chr[chr].c_str()) != 0){
			if(fread(&lchrom, sizeof(int), 1, fp_target) != 1) break;
			fread(chrom_target, sizeof(char), lchrom, fp_target);
			chrom_target[lchrom] = '\0';
			fread(&bp_target, sizeof(int), 1, fp_target);
			fread(&ancestral_target, sizeof(char), 1, fp_target);
			fread(&derived_target, sizeof(char), 1, fp_target);
			fread(&AAF_target, sizeof(int), 1, fp_target);
			fread(&DAF_target, sizeof(int), 1, fp_target);
		}

		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0){

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

					bool use = true;
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;

					if(use){
						DAF_ref    = 0;
						AAF_ref    = 0;

						while( strcmp(chrom_target, name_chr[chr].c_str()) == 0 && bp_target < bp_mut){
							if(fread(&lchrom, sizeof(int), 1, fp_target) != 1) break;
							fread(chrom_target, sizeof(char), lchrom, fp_target);
							chrom_target[lchrom] = '\0';
							fread(&bp_target, sizeof(int), 1, fp_target);
							fread(&ancestral_target, sizeof(char), 1, fp_target);
							fread(&derived_target, sizeof(char), 1, fp_target);
							fread(&AAF_target, sizeof(int), 1, fp_target);
							fread(&DAF_target, sizeof(int), 1, fp_target);
						}

						if( strcmp(chrom_target, name_chr[chr].c_str()) != 0 || bp_target != bp_mut || ancestral_target != ancestral[0] || derived_target != derived[0]){
							use = false;
						}else{
						  assert(bp_target == bp_mut);
						  os << name_chr[chr] << " "  << bp_target << " " << ancestral_target << " " << derived_target << " " << (*it_mut).age_begin << " " << (*it_mut).age_end << " " << AAF_target << " " << DAF_target << std::endl;
							//std::cerr << name_chr[chr] << " "  << bp_target << " " << ancestral_target << " " << derived_target << " " << (*it_mut).age_begin << " " << (*it_mut).age_end << " " << AAF_target << " " << DAF_target << std::endl;
						}
					}
				}
			}
		}

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
compare_tmp(cxxopts::Options& options){

	//Program options

	bool help = false;
	if(!options.count("target_tmp") || !options.count("mut") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: target_tmp, mut, output. Optional: chr." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Compare tmp." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating Colate tmp input file for " << options["target_tmp"].as<std::string>() << ".." << std::endl;

	std::uniform_real_distribution<double> dist_unif(0,1);
	std::mt19937 rng;
	int seed = std::time(0) + getpid();
	if(options.count("seed") > 0){
		seed = options["seed"].as<int>();
	}
	rng.seed(seed);

	std::string params = "20,30,10";
	if(options.count("filters")){
		params = options["filters"].as<std::string>();
	}

	std::string line;
	////////////////////////////////////////

	std::vector<std::string> name_chr, filename_mut;

	if(options.count("chr") > 0){

		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		while(getline(is_chr, line)){
			name_chr.push_back(line);
			filename_mut.push_back(options["mut"].as<std::string>() + "_chr" + line + ".mut");
		}
		is_chr.close();

	}else{

		name_chr.push_back("");
		filename_mut.push_back(options["mut"].as<std::string>());

	}

	int N_ref, N_target, DAF, bp_ref = 0, bp_target = 0, bp_mut = 0, num_used_snps = 0, num_used_snps2 = 0;
	std::string ancestral, derived;
	char ancestral_ref, derived_ref, ancestral_target, derived_target;
	Muts::iterator it_mut;

	char chrom_target[1024], chrom_ref[1024];
	int DAF_target, AAF_target, DAF_ref, AAF_ref, lchrom;

	FILE* fp_target = fopen(options["target_tmp"].as<std::string>().c_str(), "rb");
  FILE* fp_ref    = fopen(options["ref_tmp"].as<std::string>().c_str(), "rb");

  int binsize = 10e6;
	double num_mismatch = 0, num_snps = 0;
	int current_bp = 0;

	std::ofstream os(options["output"].as<std::string>());
	for(int chr = 0; chr < name_chr.size(); chr++){

		std::cerr << "parsing CHR: " << chr+1 << " / " << name_chr.size() << std::endl;

		Mutations mut;
		mut.Read(filename_mut[chr]);

		num_mismatch = 0;
		num_snps = 0;
		current_bp = mut.info[0].pos;

		while( strcmp(chrom_target, name_chr[chr].c_str()) != 0){
			if(fread(&lchrom, sizeof(int), 1, fp_target) != 1) break;
			fread(chrom_target, sizeof(char), lchrom, fp_target);
			chrom_target[lchrom] = '\0';
			fread(&bp_target, sizeof(int), 1, fp_target);
			fread(&ancestral_target, sizeof(char), 1, fp_target);
			fread(&derived_target, sizeof(char), 1, fp_target);
			fread(&AAF_target, sizeof(int), 1, fp_target);
			fread(&DAF_target, sizeof(int), 1, fp_target);
		}

		while( strcmp(chrom_ref, name_chr[chr].c_str()) != 0){
			if(fread(&lchrom, sizeof(int), 1, fp_ref) != 1) break;
			fread(chrom_ref, sizeof(char), lchrom, fp_ref);
			chrom_ref[lchrom] = '\0';
			fread(&bp_ref, sizeof(int), 1, fp_ref);
			fread(&ancestral_ref, sizeof(char), 1, fp_ref);
			fread(&derived_ref, sizeof(char), 1, fp_ref);
			fread(&AAF_ref, sizeof(int), 1, fp_ref);
			fread(&DAF_ref, sizeof(int), 1, fp_ref);
		}

		for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){

			while((*it_mut).pos > current_bp + binsize){
        os << name_chr[chr] << " " << current_bp << " " << num_mismatch << " " << num_snps << "\n";
				num_mismatch = 0;
				num_snps     = 0;
				current_bp += binsize;
			}

			if((*it_mut).flipped == 0 && (*it_mut).branch.size() == 1 && (*it_mut).age_begin < (*it_mut).age_end && (*it_mut).freq.size() >= 0){

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

					bool use = true;
					if(ancestral != "A" && ancestral != "C" && ancestral != "G" && ancestral != "T" && ancestral != "0") use = false;
					if(derived != "A" && derived != "C" && derived != "G" && derived != "T" && derived != "1") use = false;

					if(use){
						AAF_target = 0;
						DAF_target = 0;

						while( strcmp(chrom_target, name_chr[chr].c_str()) == 0 && bp_target < bp_mut){
							if(fread(&lchrom, sizeof(int), 1, fp_target) != 1) break;
							fread(chrom_target, sizeof(char), lchrom, fp_target);
							chrom_target[lchrom] = '\0';
							fread(&bp_target, sizeof(int), 1, fp_target);
							fread(&ancestral_target, sizeof(char), 1, fp_target);
							fread(&derived_target, sizeof(char), 1, fp_target);
							fread(&AAF_target, sizeof(int), 1, fp_target);
							fread(&DAF_target, sizeof(int), 1, fp_target);
						}

						DAF_ref    = 0;
						AAF_ref    = 0;
						while( strcmp(chrom_ref, name_chr[chr].c_str()) == 0 && bp_ref < bp_mut){
							if(fread(&lchrom, sizeof(int), 1, fp_ref) != 1) break;
							fread(chrom_ref, sizeof(char), lchrom, fp_ref);
							chrom_ref[lchrom] = '\0';
							fread(&bp_ref, sizeof(int), 1, fp_ref);
							fread(&ancestral_ref, sizeof(char), 1, fp_ref);
							fread(&derived_ref, sizeof(char), 1, fp_ref);
							fread(&AAF_ref, sizeof(int), 1, fp_ref);
							fread(&DAF_ref, sizeof(int), 1, fp_ref);
						}

						if( strcmp(chrom_target, name_chr[chr].c_str()) != 0 || bp_target != bp_mut || ancestral_target != ancestral[0] || derived_target != derived[0]){
							use = false;
						}else if(strcmp(chrom_ref, name_chr[chr].c_str()) != 0 || bp_ref != bp_mut || ancestral_ref != ancestral[0] || derived_ref != derived[0]){
							use = false;
						}else{
							assert(bp_target == bp_mut);
              assert(bp_ref    == bp_mut);
							//now sample a read/genotype from each and add to pairwise mismatch
           
							if(0){
								int target_sampled = 0;
								int ref_sampled    = 0;

								if(dist_unif(rng) < DAF_target/(AAF_target + DAF_target)){
									target_sampled = 1;
								}
								if(dist_unif(rng) < DAF_ref/(AAF_ref + DAF_ref)){
									ref_sampled = 1;
								}
								num_mismatch += fabs(target_sampled - ref_sampled);
							}

							num_mismatch += (DAF_target * (N_ref - DAF_ref) + (N_target - DAF_target) * DAF_ref)/((double) N_target * N_ref);
							num_snps++;

						}
					}
				}
			}
		}

		//print also at end of chromosome
		os << name_chr[chr] << " " << current_bp << " " << num_mismatch << " " << num_snps << "\n";
		num_mismatch = 0;
		num_snps     = 0;
		current_bp += binsize;

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




//////////////

void
GetConditionalCoalescentRate(Tree& tree, Sample& sample, std::vector<Leaves>& leaves, float factor, std::vector<int>& focal_haps, std::vector<int>& conditional_haps, std::vector<float>& epoch, std::vector<float>& epoch_focal, std::vector<CollapsedMatrix<float>>& coalescence_rate_num, std::vector<CollapsedMatrix<float>>& coalescence_rate_denom){

	for(std::vector<int>::iterator it_focal = focal_haps.begin(); it_focal != focal_haps.end(); it_focal++){
		for(std::vector<int>::iterator it_cond = conditional_haps.begin(); it_cond != conditional_haps.end(); it_cond++){

			Node node = tree.nodes[*it_focal];
			int parent = (*tree.nodes[*it_focal].parent).label;
			int child;
			float coord = 0.0;
			int ep_start = 0, ep = 0, ep_init = 0;
			float coal_age  = 0.0;
			float lower_age = 0.0;

			bool use = false;
			if(*it_cond == -1) use = true;
			while(parent < tree.nodes.size()){

				if(!use){
					for(std::vector<int>::iterator it_mem = leaves[node.label].member.begin(); it_mem != leaves[node.label].member.end(); it_mem++){
						if(*it_cond == *it_mem){
							coal_age = coord;
							use = true;
							break;
						}
					}
					ep_start = 0;
					if(epoch_focal[ep_start] < coal_age){
						while(epoch_focal[ep_start] < coal_age){
							ep_start++;
							if(ep_start == epoch_focal.size()) break;
						}
						ep_start--;
					}
					assert(epoch_focal[ep_start] <= coal_age);
					assert(epoch_focal[ep_start + 1] > coal_age);
				}

				coord += node.branch_length;   
				if(use){

					ep_init = 0;
					if(coal_age > epoch[ep_init]){
						while(coal_age > epoch[ep_init]){
							ep_init++; 
						}
						ep_init--;
					}

					child = (*tree.nodes[parent].child_left).label; 
					if(child == node.label) child = (*tree.nodes[parent].child_right).label;

					for(std::vector<int>::iterator it_other = leaves[child].member.begin(); it_other != leaves[child].member.end(); it_other++){

						lower_age = coal_age;
						ep = ep_init;
						while(coord > epoch[ep+1]){
							coalescence_rate_denom[ep_start][ep][sample.group_of_haplotype[*it_other]] += factor * (epoch[ep+1] - lower_age);
							ep++;
							if(ep == epoch.size()) break;
							lower_age = epoch[ep]; 
						}
						assert(coord < epoch[ep+1]);
						coalescence_rate_denom[ep_start][ep][sample.group_of_haplotype[*it_other]] += factor * (coord - lower_age);
						coalescence_rate_num[ep_start][ep][sample.group_of_haplotype[*it_other]]   += factor;

					}
				}

				node = (*node.parent);
				if(node.label == tree.nodes.size()-1) break;
				parent = (*tree.nodes[node.label].parent).label; 
			} 

		}
	}

}

void
GetConditionalCoalescentRate(Tree& tree, Sample& sample, std::vector<Leaves>& leaves, float factor, std::vector<int>& focal_haps, std::vector<int>& conditional_haps, std::vector<float>& epoch, std::vector<float>& epoch_focal, std::vector<double>& sample_ages, std::vector<CollapsedMatrix<float>>& coalescence_rate_num, std::vector<CollapsedMatrix<float>>& coalescence_rate_denom){


	for(std::vector<int>::iterator it_focal = focal_haps.begin(); it_focal != focal_haps.end(); it_focal++){
		for(std::vector<int>::iterator it_cond = conditional_haps.begin(); it_cond != conditional_haps.end(); it_cond++){

			Node node = tree.nodes[*it_focal];
			int parent = (*tree.nodes[*it_focal].parent).label, child;
			double age = sample_ages[*it_focal];
			float coal_age = age;
			float coord = age;

			int ep_coal = 0, ep_start = 0, ep = 0, ep_init = 0;

			float lower_age = age;
			bool use = false;
			if(*it_cond == -1) use = true;
			while(parent < tree.nodes.size()){

				if(!use){
					for(std::vector<int>::iterator it_mem = leaves[node.label].member.begin(); it_mem != leaves[node.label].member.end(); it_mem++){
						if(*it_cond == *it_mem){
							coal_age = coord;
							use = true;
							break;
						}
					}
					ep_coal = 0;
					if(epoch_focal[ep_coal] < coord){
						while(epoch_focal[ep_coal] < coord){
							ep_coal++;
							if(ep_coal == epoch_focal.size()) break;
						}
						ep_coal--;
					}
				}
				coord += node.branch_length;

				if(use){
					child = (*tree.nodes[parent].child_left).label; 
					if(child == node.label) child = (*tree.nodes[parent].child_right).label;

					ep_init = 0;
					if(coal_age > epoch[ep_init]){
						while(coal_age > epoch[ep_init]){
							ep_init++; 
						}
						ep_init--;
					}

					for(std::vector<int>::iterator it_other = leaves[child].member.begin(); it_other != leaves[child].member.end(); it_other++){

						lower_age = std::max(age, sample_ages[*it_other]);
						lower_age = std::max(lower_age, coal_age);

						if(lower_age <= coal_age){
							ep_start = ep_coal;
							ep = ep_init;
						}else{
							ep_start = ep_coal;
							if(epoch_focal[ep_start] < lower_age){
								while(epoch_focal[ep_start] < lower_age){
									ep_start++;
									if(ep_start == epoch_focal.size()) break;
								}
								ep_start--;
							}
							ep = ep_init;
							if(epoch[ep] < lower_age){
								while(epoch[ep] < lower_age){
									ep++;
									if(ep == epoch.size()) break;
								}
								ep--;
							}
						}

						while(coord > epoch[ep+1]){
							coalescence_rate_denom[ep_start][ep][sample.group_of_haplotype[*it_other]] += factor * (epoch[ep+1] - lower_age);
							ep++;
							if(ep == epoch.size()) break;
							lower_age = epoch[ep]; 
						}
						coalescence_rate_denom[ep_start][ep][sample.group_of_haplotype[*it_other]] += factor * (coord - lower_age);
						coalescence_rate_num[ep_start][ep][sample.group_of_haplotype[*it_other]]   += factor;

					}
				}


				node = (*node.parent);
				if(node.label == tree.nodes.size()-1) break;
				parent = (*tree.nodes[node.label].parent).label; 
			} 

		}
	}

}

int
CondCoalRates(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("input") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output. Optional: years_per_gen, dist, bins." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Reads .anc file and calculates conditional pairwise coalescence rate. Output is bin file. Use SummarizeCoalesecntRate to obtain coalescence rates." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating conditional coalescence rate for " << options["input"].as<std::string>() << " ..." << std::endl;

	////////////////////////
	//read in anc file

	std::string line;
	std::vector<std::string> filename_chr;
	if(options.count("chr")){
		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		while(getline(is_chr, line)){
			filename_chr.push_back(line);
		}
		is_chr.close();
	}else{
		filename_chr.push_back("NA");
	}

	AncMutIterators ancmut1;
	if(options.count("dist")){
		if(filename_chr[0] == "NA"){
			ancmut1.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut", options["dist"].as<std::string>());
		}else{
			ancmut1.OpenFiles(options["input"].as<std::string>() + "_chr" + filename_chr[0] + ".anc", options["input"].as<std::string>() + "_chr" + filename_chr[0] + ".mut", options["dist"].as<std::string>() + "_chr" + filename_chr[0] + ".dist");
		} 
	}else{  
		if(filename_chr[0] == "NA"){
			ancmut1.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
		}else{
			ancmut1.OpenFiles(options["input"].as<std::string>() + "_chr" + filename_chr[0] + ".anc", options["input"].as<std::string>() + "_chr" + filename_chr[0] + ".mut");
		} 
	} 

	int N = ancmut1.NumTips();
	int L = ancmut1.NumSnps();
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0; 
	Data data(N,L); 

	//calculate epoch times

	float years_per_gen = 28.0;
	if(options.count("years_per_gen")){
		years_per_gen = options["years_per_gen"].as<float>();
	}

	int num_epochs;
	std::vector<float> epochs;
	float log_10 = std::log(10);
	if(options.count("bins")){

		double log_age = std::log(0);
		double age = 0;

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
				if(ep == 1 && age != 0.0) epochs.push_back(age);
				epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
				ep++;
			}
			epoch_boundary += epoch_step;
		}
		epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
		epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();

	}else{

		num_epochs = 31;
		epochs.resize(num_epochs);
		epochs[0] = 0.0;
		epochs[1] = 1e3/years_per_gen;
		for(int e = 2; e < num_epochs-1; e++){
			epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
		}
		epochs[num_epochs-1] = 1e8/years_per_gen;

	}

	std::vector<float> epochs_focal = {0,1e4,5e4,1e5,1e6,1e7};
	for(std::vector<float>::iterator it_ep = epochs_focal.begin(); it_ep != epochs_focal.end(); it_ep++){
		*it_ep /= years_per_gen;
	}

	std::mt19937 rng;
	int seed = std::time(0) + getpid();
	if(options.count("seed") > 0){
		seed = options["seed"].as<int>();
	}
	rng.seed(seed);

	Sample sample;
	sample.Read(options["poplabels"].as<std::string>());

	std::string groups = options["groups"].as<std::string>();
	//format focal_group;conditional_group
	std::string g1,g2;
	int i = 0;
	while(groups[i] != ','){
		g1 += groups[i];
		i++;
		if(i == groups.size()) break;
	}
	i++;
	while(i < groups.size()){
		g2 += groups[i];
		i++;
	}

	std::cerr << g1 << "|" << g2 << std::endl;

	std::vector<int> focal_haps;
	std::vector<int> conditional_haps;
	for(int i = 0; i < data.N; i++){
		if(sample.groups[sample.group_of_haplotype[i]] == g1) focal_haps.push_back(i);
		if(sample.groups[sample.group_of_haplotype[i]] == g2) conditional_haps.push_back(i);
	} 
	if(conditional_haps.size() == 0) conditional_haps.push_back(-1);
	if(focal_haps.size() == 0){
		std::cerr << "Error: groups not found" << std::endl;
		exit(1);
	}


	int bin_size = 30e6;
	int num_blocks = (int)((ancmut1.NumPos())/bin_size)+2;

	std::vector<std::vector<CollapsedMatrix<float>>> coalescence_rate_num(num_blocks), coalescence_rate_denom(num_blocks);
	for(int bin = 0; bin < num_blocks; bin++){
		coalescence_rate_num[bin].resize(num_epochs);
		coalescence_rate_denom[bin].resize(num_epochs);
		for(int e = 0; e < epochs_focal.size(); e++){
			coalescence_rate_num[bin][e].resize(num_epochs, sample.groups.size());
			std::fill(coalescence_rate_num[bin][e].vbegin(), coalescence_rate_num[bin][e].vend(), 0.0);
			coalescence_rate_denom[bin][e].resize(num_epochs, sample.groups.size());
			std::fill(coalescence_rate_denom[bin][e].vbegin(), coalescence_rate_denom[bin][e].vend(), 0.0);
		}
	}

	int bin = 0, chr_bin = 0;
	for(int chr = 0; chr < filename_chr.size(); chr++){

		std::cerr << "CHR: " << filename_chr[chr] << std::endl;


		AncMutIterators ancmut;
		if(options.count("dist")){
			if(filename_chr[0] == "NA"){
				ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut", options["dist"].as<std::string>());
			}else{
				ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + filename_chr[chr] + ".anc", options["input"].as<std::string>() + "_chr" + filename_chr[chr] + ".mut", options["dist"].as<std::string>() + "_chr" + filename_chr[chr] + ".dist");
			} 
		}else{  
			if(filename_chr[0] == "NA"){
				ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
			}else{
				ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + filename_chr[chr] + ".anc", options["input"].as<std::string>() + "_chr" + filename_chr[chr] + ".mut");
			} 
		} 

		fasta mask;
		double cutoff = 0.9;
		if(options.count("mask")){
			if(filename_chr[0] == "NA"){
				mask.Read(options["mask"].as<std::string>());
			}else{
				mask.Read(options["mask"].as<std::string>() + "_chr" + filename_chr[chr] + ".fa");
			}
		}


		bin = chr_bin;

		int extra_blocks = (int)(ancmut.NumPos()/bin_size)+2; 

		if(bin + extra_blocks > num_blocks){
			int bin_tmp = num_blocks;
			num_blocks += extra_blocks;
			coalescence_rate_num.resize(num_blocks);
			coalescence_rate_denom.resize(num_blocks);
			for(; bin_tmp < num_blocks; bin_tmp++){
				coalescence_rate_num[bin_tmp].resize(num_epochs);
				coalescence_rate_denom[bin_tmp].resize(num_epochs);
				for(int e = 0; e < epochs_focal.size(); e++){
					coalescence_rate_num[bin_tmp][e].resize(num_epochs, sample.groups.size());
					std::fill(coalescence_rate_num[bin_tmp][e].vbegin(), coalescence_rate_num[bin_tmp][e].vend(), 0.0);
					coalescence_rate_denom[bin_tmp][e].resize(num_epochs, sample.groups.size());
					std::fill(coalescence_rate_denom[bin_tmp][e].vbegin(), coalescence_rate_denom[bin_tmp][e].vend(), 0.0);
				}
			}
		}

		////////////////////////////////
		//Pairwise coalescence rate
		//In each tree, find the coalescent time. Then update count_per_epoch and coalescent_time_in_epoch. 


		if(ancmut.sample_ages.size() > 0){
			float factor = 0.0;
			num_bases_tree_persists = 0.0;
			double num_passing = 1.0;
			while(num_bases_tree_persists >= 0.0){
				num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
				bin = (int)((*it_mut).pos/bin_size) + chr_bin;

				num_passing = 1.0;
				if(options.count("mask") > 0){

					int tree_index = (*it_mut).tree;
					int pos_start = (*it_mut).pos;
					int pos_end = pos_start;
					if(it_mut != ancmut.mut_end()){
						while((*it_mut).tree == tree_index){
							pos_end = (*it_mut).pos;
							it_mut++;
							if(it_mut == ancmut.mut_end()) break;
						}
					}
					num_passing = 0.0;
					if(pos_start < mask.seq.size() && pos_end < mask.seq.size()){
						for(int bp = pos_start; bp < pos_end; bp++){
							if(mask.seq[bp-1] == 'P') num_passing++; 
						}
					}

					if(pos_end - pos_start + 1 <= 0){ 
						num_passing = 0.0;
					}else{
						num_passing /= (pos_end - pos_start + 1);
					}

				}

				if(num_passing >= cutoff){
					std::vector<Leaves> leaves;
					mtr.tree.FindAllLeaves(leaves);
					factor = num_bases_tree_persists;
					GetConditionalCoalescentRate(mtr.tree, sample, leaves, factor, focal_haps, conditional_haps, epochs, epochs_focal, ancmut.sample_ages, coalescence_rate_num[bin], coalescence_rate_denom[bin]);
				}

			}  

		}else{
			float factor = 0.0;
			num_bases_tree_persists = 0.0;
			double num_passing = 1.0;
			while(num_bases_tree_persists >= 0.0){
				num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
				bin = (int)((*it_mut).pos/bin_size) + chr_bin;

				num_passing = 1.0;
				num_passing = 1.0;
				if(options.count("mask") > 0){

					int tree_index = (*it_mut).tree;
					int pos_start = (*it_mut).pos;
					int pos_end = pos_start;
					if(it_mut != ancmut.mut_end()){
						while((*it_mut).tree == tree_index){
							pos_end = (*it_mut).pos;
							it_mut++;
							if(it_mut == ancmut.mut_end()) break;
						}
					}
					num_passing = 0.0;
					if(pos_start < mask.seq.size() && pos_end < mask.seq.size()){
						for(int bp = pos_start; bp < pos_end; bp++){
							if(mask.seq[bp-1] == 'P') num_passing++; 
						}
					}

					if(pos_end - pos_start + 1 <= 0){ 
						num_passing = 0.0;
					}else{
						num_passing /= (pos_end - pos_start + 1);
					}

				}

				if(num_passing >= cutoff){
					std::vector<Leaves> leaves;
					mtr.tree.FindAllLeaves(leaves);
					factor = num_bases_tree_persists;
					GetConditionalCoalescentRate(mtr.tree, sample, leaves, factor, focal_haps, conditional_haps, epochs, epochs_focal, coalescence_rate_num[bin], coalescence_rate_denom[bin]);
				}
			}  
		}

		chr_bin = bin+1;

	}

	num_blocks = bin+1;
	coalescence_rate_num.resize(num_blocks);
	coalescence_rate_denom.resize(num_blocks);

	int num_bootstrap = 1;
	if(options.count("num_bootstraps") > 0) num_bootstrap = options["num_bootstraps"].as<int>();

	std::uniform_int_distribution<int> dist_blocks(0,num_blocks-1);
	std::vector<int> blocks(num_blocks);


	//output  
	std::ofstream os(options["output"].as<std::string>());
	os << "boot lineage_epoch epoch.start group rate" << std::endl;
	for(int iter = 0; iter < num_bootstrap; iter++){

		if(num_bootstrap == 1){
			std::fill(blocks.begin(), blocks.end(), 1.0);
		}else{
			std::fill(blocks.begin(), blocks.end(), 0.0);
			for(int block = 0; block < num_blocks; block++){
				blocks[dist_blocks(rng)] += 1.0;
			}
		}

		std::vector<CollapsedMatrix<float>> res_num(num_epochs), res_denom(num_epochs);
		for(int e = 0; e < epochs_focal.size(); e++){
			res_num[e].resize(num_epochs, sample.groups.size());
			std::fill(res_num[e].vbegin(), res_num[e].vend(), 0.0);
			res_denom[e].resize(num_epochs, sample.groups.size());
			std::fill(res_denom[e].vbegin(), res_denom[e].vend(), 0.0);
		}

		for(int block = 0; block < num_blocks; block++){

			for(int e = 0; e < epochs_focal.size(); e++){
				std::vector<float>::iterator it_coal_num = coalescence_rate_num[block][e].vbegin();
				for(std::vector<float>::iterator it_num = res_num[e].vbegin(); it_num != res_num[e].vend(); it_num++){
					(*it_num) += blocks[block] * (*it_coal_num);
					it_coal_num++;
				}
				std::vector<float>::iterator it_coal_denom = coalescence_rate_denom[block][e].vbegin();
				for(std::vector<float>::iterator it_denom = res_denom[e].vbegin(); it_denom != res_denom[e].vend(); it_denom++){
					(*it_denom) += blocks[block] * (*it_coal_denom);
					it_coal_denom++;
				}
			}

		}

		for(int ep1 = 0; ep1 < epochs_focal.size(); ep1++){
			for(int ep2 = 0; ep2 < num_epochs; ep2++){
				for(int i = 0; i < sample.groups.size(); i++){
					os << iter << " " << epochs_focal[ep1] << " " << epochs[ep2] << " " << sample.groups[i] << " " << res_num[ep1][ep2][i]/res_denom[ep1][ep2][i] << std::endl;
				}
			}
		}

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

	return 0;

}



