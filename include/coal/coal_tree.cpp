#include "coal_tree.hpp"

coal_tree::coal_tree(std::vector<double> epochs, int num_bootstrap, int block_size, AncMutIterators& ancmut): epochs(epochs), num_bootstrap(num_bootstrap), block_size(block_size){

	N = ancmut.NumTips();
	N_total = 2*N-1;
	coords.resize(N_total);
	num_lins.resize(N_total);
	sorted_indices.resize(N_total);

	num_trees = ancmut.NumTrees();
	num_blocks = num_trees/((double) block_size) + 1;
	num_epochs = epochs.size();

	num.resize(num_blocks);
	denom.resize(num_blocks);
	for(it1_num = num.begin(); it1_num != num.end(); it1_num++){
		(*it1_num).resize(num_epochs);
		std::fill((*it1_num).begin(), (*it1_num).end(), 0.0);
	}
	for(it1_denom = denom.begin(); it1_denom != denom.end(); it1_denom++){
		(*it1_denom).resize(num_epochs);
		std::fill((*it1_denom).begin(), (*it1_denom).end(), 0.0);
	}

	rng.seed(1);
	std::uniform_int_distribution<int> d(0, num_blocks);
	std::vector<int> tmp(num_blocks);
	blocks.resize(num_bootstrap);
	for(it1_blocks = blocks.begin(); it1_blocks != blocks.end(); it1_blocks++){

		(*it1_blocks).resize(num_blocks);
		for(it2_blocks = (*it1_blocks).begin(); it2_blocks != (*it1_blocks).end(); it2_blocks++){
			*it2_blocks = d(rng);
		}
		std::sort((*it1_blocks).begin(), (*it1_blocks).end());
		it2_blocks = (*it1_blocks).begin();
		for(int i = 0; i < num_blocks; i++){	
			tmp[i] = 0;
			if(it2_blocks != (*it1_blocks).end()){
				if(*it2_blocks == i){
					while(*it2_blocks == i){
						tmp[i]++;
						it2_blocks++;
						if(it2_blocks == (*it1_blocks).end()) break;
					}
				}
			}
			//std::cerr << tmp[i] << " ";
		}
		//std::cerr << std::endl;
		*it1_blocks = tmp;

	}
	num_boot.resize(num_bootstrap);
	denom_boot.resize(num_bootstrap);
	it1_num = num_boot.begin();
	it1_denom = denom_boot.begin();
	for(;it1_num != num_boot.end();){
    (*it1_num).resize(num_epochs);
		(*it1_denom).resize(num_epochs);
		it1_num++;
		it1_denom++;
	}

	it1_num   = num.begin();
	it1_denom = denom.begin();
	current_block = 0;
	count_trees   = 0;

}

void 
coal_tree::populate(Tree& tree){

	tree.GetCoordinates(coords);

	std::size_t m1(0);
	std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
	std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
			return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

	int lins = 0;
	float age = coords[*sorted_indices.begin()];
	it_sorted_indices      = sorted_indices.begin();
	it_sorted_indices_prev = it_sorted_indices;
	it_num_lins            = num_lins.begin();
	for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
		if(coords[*it_sorted_indices] > age){
			while(coords[*it_sorted_indices_prev] == age){
				*it_num_lins = lins;
				it_num_lins++;
				it_sorted_indices_prev++;
			}
			age = coords[*it_sorted_indices_prev];
			assert(it_sorted_indices_prev == it_sorted_indices);
		}
		if(*it_sorted_indices < N){
			lins++;
		}else{
			lins--;
		}
		assert(lins >= 1);
	}
	while(coords[*it_sorted_indices_prev] == age){
		*it_num_lins = lins;
		it_num_lins++;
		it_sorted_indices_prev++;
		if(it_num_lins == num_lins.end()) break;
	}
	//populate num and denom
	if(count_trees == block_size){
		current_block++;
		count_trees = 0;
		it1_num++;
		it1_denom++;
		assert(current_block < num_blocks);
	}
	std::sort(coords.begin(), coords.end());

	/*
	for(int i = 0; i < 2*N-1; i++){
		std::cerr << i << " " << coords[i] << " " << num_lins[i] << std::endl;
	}
	std::cerr << std::endl;
	*/

	it2_num           = (*it1_num).begin();
	it2_denom         = (*it1_denom).begin();
	
	it_num_lins       = num_lins.begin();
	it_coords_next    = std::next(coords.begin(), 1);
	it_sorted_indices = std::next(sorted_indices.begin(),1);
	it_epochs         = std::next(epochs.begin(),1);
	double current_lower_age = epochs[0];
	for(;it_epochs != epochs.end();){

		while(*it_coords_next <= *it_epochs){
			if(*it_sorted_indices >= N) (*it2_num)++;
			*it2_denom += (*it_num_lins) * (*it_num_lins-1)/2.0 * (*it_coords_next - current_lower_age)/num_trees;
			current_lower_age = *it_coords_next;
			it_num_lins++;
			it_coords_next++;
			it_sorted_indices++;
			if(it_coords_next == coords.end()) break;
		}
		if(it_coords_next == coords.end()) break;
		*it2_denom += (*it_num_lins) * (*it_num_lins-1)/2.0 * (*it_epochs - current_lower_age)/num_trees;
		current_lower_age = *it_epochs;
		it_epochs++;
		it2_num++;
		it2_denom++;

	}
  assert(*it_num_lins == 1);
	count_trees++;

}

void 
coal_tree::Dump(std::string filename){

	//populate num_boot and num_denom (vector of size num_bootstrap)
	//using num, denom, blocks

	std::vector<std::vector<double>>::iterator it1_num_boot = num_boot.begin(), it1_denom_boot = denom_boot.begin();
	std::vector<double>::iterator it2_num_boot, it2_denom_boot;
	for(it1_blocks = blocks.begin(); it1_blocks != blocks.end(); it1_blocks++){

		it2_num_boot   = (*it1_num_boot).begin();
		it2_denom_boot = (*it1_denom_boot).begin();
		for(;it2_num_boot != (*it1_num_boot).end();){
			*it2_num_boot   = 0.0;
			*it2_denom_boot = 0.0;
			it2_num_boot++;
			it2_denom_boot++;
		}

		it1_num   = num.begin();
		it1_denom = denom.begin();
		assert(num.size() == (*it1_blocks).size());
		for(it2_blocks = (*it1_blocks).begin(); it2_blocks != (*it1_blocks).end(); it2_blocks++){

			if(*it2_blocks > 0){

				assert((*it1_num).size() == (*it1_num_boot).size());
				it2_num        = (*it1_num).begin();
				it2_denom      = (*it1_denom).begin();
				it2_num_boot   = (*it1_num_boot).begin();
				it2_denom_boot = (*it1_denom_boot).begin();
				for(;it2_num != (*it1_num).end();){
					//std::cerr << *it2_blocks << " " << *it2_num << " " << *it2_denom << std::endl;
					*it2_num_boot   += *it2_blocks * *it2_num;
					*it2_denom_boot += *it2_blocks * *it2_denom;
					it2_num++;
					it2_denom++;
					it2_num_boot++;
					it2_denom_boot++;
				}

			}
			//std::cerr << "done" << std::endl;

			it1_num++;
			it1_denom++;
		}

		it1_num_boot++;
		it1_denom_boot++;
	}

	std::ofstream os(filename);
	for(int i = 0; i < num_bootstrap; i++){
    os << i << " ";
	}
	os << "\n";
	for(it_epochs = epochs.begin(); it_epochs != epochs.end(); it_epochs++){
		os << *it_epochs << " ";
	}
	os << "\n";
	it1_num_boot = num_boot.begin();
	it1_denom_boot = denom_boot.begin();
	int i = 0;
	for(; it1_num_boot != num_boot.end();){

		os << "0 " << i << " ";
		it2_num_boot = (*it1_num_boot).begin();
		it2_denom_boot = (*it1_denom_boot).begin();
		for(;it2_num_boot != (*it1_num_boot).end();){
			os << *it2_num_boot/((*it2_denom_boot)*num_trees) << " ";
			it2_num_boot++;
			it2_denom_boot++;
		}
		os << "\n";

		i++;
		it1_num_boot++;
		it1_denom_boot++;
	}
	os.close();

}
