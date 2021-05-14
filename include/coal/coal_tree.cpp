#include "coal_tree.hpp"

////////////////////////////

coal_tree::coal_tree(std::vector<double> epochs, int num_bootstrap, int block_size): epochs(epochs), num_bootstrap(num_bootstrap), block_size(block_size){

  num_blocks = 0;
	num_epochs = epochs.size();

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
coal_tree::update_ancmut(AncMutIterators& ancmut){

	N = ancmut.NumTips();
	N_total = 2*N-1;
	coords.resize(N_total);
	num_lins.resize(N_total);
	sorted_indices.resize(N_total);

	num_trees = ancmut.NumTrees();
  int num_blocks_prev = num_blocks;
  num_blocks += num_trees/((double) block_size) + 1;

	num.resize(num_blocks);
	denom.resize(num_blocks);
	for(it1_num = std::next(num.begin(),num_blocks_prev); it1_num != num.end(); it1_num++){
		(*it1_num).resize(num_epochs);
		std::fill((*it1_num).begin(), (*it1_num).end(), 0.0);
	}
	for(it1_denom = std::next(denom.begin(),num_blocks_prev); it1_denom != denom.end(); it1_denom++){
		(*it1_denom).resize(num_epochs);
		std::fill((*it1_denom).begin(), (*it1_denom).end(), 0.0);
	}

	it1_num   = std::next(num.begin(), num_blocks_prev);
	it1_denom = std::next(denom.begin(), num_blocks_prev);
  current_block = num_blocks_prev;
  count_trees = 0;

}

void 
coal_tree::populate(Tree& tree, double num_bases_tree_persists){

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
			if(*it_sorted_indices >= N) (*it2_num) += num_bases_tree_persists/1e9;
			*it2_denom += num_bases_tree_persists * (*it_num_lins) * (*it_num_lins-1)/2.0 * (*it_coords_next - current_lower_age)/1e9;
			current_lower_age = *it_coords_next;
			it_num_lins++;
			it_coords_next++;
			it_sorted_indices++;
			if(it_coords_next == coords.end()) break;
		}
		if(it_coords_next == coords.end()) break;
		*it2_denom += num_bases_tree_persists * (*it_num_lins) * (*it_num_lins-1)/2.0 * (*it_epochs - current_lower_age)/1e9;
		current_lower_age = *it_epochs;
		it_epochs++;
		it2_num++;
		it2_denom++;

	}
  assert(*it_num_lins == 1);
	count_trees++;

}

void
coal_tree::init_bootstrap(){

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

}

void 
coal_tree::Dump(std::string filename){

	//populate num_boot and num_denom (vector of size num_bootstrap)
	//using num, denom, blocks

  init_bootstrap();

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
			std::cerr << *it2_num_boot << " " << (*it2_denom_boot) << std::endl;
			os << *it2_num_boot/((*it2_denom_boot)) << " ";
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




////////////////////////////

coal_LA::coal_LA(std::vector<double> epochs, int num_bootstrap, int block_size, int num_groups): epochs(epochs), num_bootstrap(num_bootstrap), block_size(block_size), num_groups(num_groups){

	num_blocks = 0;
	num_epochs = epochs.size();

	num_boot.resize(num_bootstrap);
	denom_boot.resize(num_bootstrap);
	it1_num = num_boot.begin();
	it1_denom = denom_boot.begin();
	for(;it1_num != num_boot.end();){
		(*it1_num).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
      (*it1_num)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
        (*it1_num)[i][j].resize(num_epochs);
				std::fill((*it1_num)[i][j].begin(), (*it1_num)[i][j].end(), 0.0);
			}
	  }
		(*it1_denom).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
			(*it1_denom)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
				(*it1_denom)[i][j].resize(num_epochs);
				std::fill((*it1_denom)[i][j].begin(), (*it1_denom)[i][j].end(), 0.0);
			}
		}
		it1_num++;
		it1_denom++;
	}

	it1_num   = num.begin();
	it1_denom = denom.begin();
	current_block = 0;
	count_trees   = 0;

}

coal_LA::coal_LA(std::vector<double> epochs, int num_bootstrap, int block_size, int num_groups, AncMutIterators& ancmut): epochs(epochs), num_bootstrap(num_bootstrap), block_size(block_size), num_groups(num_groups){

	N = ancmut.NumTips();
	N_total = 2*N-1;
	coords.resize(N_total);
	sorted_indices.resize(N_total);

	num_trees = ancmut.NumTrees();
	num_blocks = num_trees/((double) block_size) + 1;
	num_epochs = epochs.size();

	num.resize(num_blocks);
	denom.resize(num_blocks);
	for(it1_num = num.begin(); it1_num != num.end(); it1_num++){
		(*it1_num).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
      (*it1_num)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
        (*it1_num)[i][j].resize(num_epochs);
				std::fill((*it1_num)[i][j].begin(), (*it1_num)[i][j].end(), 0.0);
			}
	  }
	}
	for(it1_denom = denom.begin(); it1_denom != denom.end(); it1_denom++){
		(*it1_denom).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
      (*it1_denom)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
        (*it1_denom)[i][j].resize(num_epochs);
				std::fill((*it1_denom)[i][j].begin(), (*it1_denom)[i][j].end(), 0.0);
			}
	  }
	}

	num_boot.resize(num_bootstrap);
	denom_boot.resize(num_bootstrap);
	it1_num = num_boot.begin();
	it1_denom = denom_boot.begin();
	for(;it1_num != num_boot.end();){
		(*it1_num).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
      (*it1_num)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
        (*it1_num)[i][j].resize(num_epochs);
				std::fill((*it1_num)[i][j].begin(), (*it1_num)[i][j].end(), 0.0);
			}
	  }
		(*it1_denom).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
			(*it1_denom)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
				(*it1_denom)[i][j].resize(num_epochs);
				std::fill((*it1_denom)[i][j].begin(), (*it1_denom)[i][j].end(), 0.0);
			}
		}
		it1_num++;
		it1_denom++;
	}

	it1_num   = num.begin();
	it1_denom = denom.begin();
	current_block = 0;
	count_trees   = 0;

}

void
coal_LA::update_ancmut(AncMutIterators& ancmut){

	N = ancmut.NumTips();
	N_total = 2*N-1;
	coords.resize(N_total);
	sorted_indices.resize(N_total);

	num_trees = ancmut.NumTrees();
	int num_blocks_prev = num_blocks;
	num_blocks += num_trees/((double) block_size) + 1;

	num.resize(num_blocks);
	denom.resize(num_blocks);
	for(it1_num = std::next(num.begin(),num_blocks_prev); it1_num != num.end(); it1_num++){
		(*it1_num).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
			(*it1_num)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
				(*it1_num)[i][j].resize(num_epochs);
				std::fill((*it1_num)[i][j].begin(), (*it1_num)[i][j].end(), 0.0);
			}
		}
	}
	for(it1_denom = std::next(denom.begin(),num_blocks_prev); it1_denom != denom.end(); it1_denom++){
		(*it1_denom).resize(num_groups);
		for(int i = 0; i < num_groups; i++){
			(*it1_denom)[i].resize(num_groups);
			for(int j = 0; j < num_groups; j++){
				(*it1_denom)[i][j].resize(num_epochs);
				std::fill((*it1_denom)[i][j].begin(), (*it1_denom)[i][j].end(), 0.0);
			}
		}
	}

	it1_num   = std::next(num.begin(), num_blocks_prev);
	it1_denom = std::next(denom.begin(), num_blocks_prev);
	current_block = num_blocks_prev;
	count_trees = 0;

}

void 
coal_LA::populate(Tree& tree, double num_bases_tree_persists, std::vector<int>& group, bool new_tree = false){

	tree.GetCoordinates(coords);
  tree.FindAllLeaves(desc);

	//populate num and denom
	if(count_trees == block_size){
		current_block++;
		count_trees = 0;
		it1_num++;
		it1_denom++;
		assert(current_block < num_blocks);
	}

	std::size_t m1(0);
	std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
	std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
			return std::tie(coords[i1],i1) < std::tie(coords[i2],i2); } );

  int ep = 1;
	it_sorted_indices      = sorted_indices.begin();
	std::vector<double> denom_tmpl(epochs.size(), 0.0);
	double lower_boundary  = 0.0;
	for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){

		while(coords[*it_sorted_indices] > epochs[ep]){
			denom_tmpl[ep-1] += (epochs[ep] - lower_boundary) * num_bases_tree_persists/1e9;
			lower_boundary = epochs[ep];
			ep++;
		}
    assert(coords[*it_sorted_indices] >= lower_boundary);
    assert(coords[*it_sorted_indices] <= epochs[ep]);

		denom_tmpl[ep-1]      += (coords[*it_sorted_indices] - lower_boundary) * num_bases_tree_persists/1e9;
		lower_boundary         = coords[*it_sorted_indices];

    if(desc[*it_sorted_indices].member.size() > 1){
      int child1 = (*tree.nodes[*it_sorted_indices].child_left).label;
      int child2 = (*tree.nodes[*it_sorted_indices].child_right).label;

      for(std::vector<int>::iterator it_mem1 = desc[child1].member.begin(); it_mem1 != desc[child1].member.end(); it_mem1++){
        for(std::vector<int>::iterator it_mem2 = desc[child2].member.begin(); it_mem2 != desc[child2].member.end(); it_mem2++){
          if(*it_mem1 != *it_mem2){
            int group1 = group[*it_mem1];
            int group2 = group[*it_mem2];
            if(group2 > group1){
              int tmp = group1;
              group1 = group2;
              group2 = tmp;
            }
            (*it1_num)[group1][group2][ep-1] += num_bases_tree_persists/1e9;
            std::vector<double>::iterator it_denom_tmpl = denom_tmpl.begin();
            for(std::vector<double>::iterator it = (*it1_denom)[group1][group2].begin(); it != (*it1_denom)[group1][group2].end(); it++){
              *it += *it_denom_tmpl;
              if(*it_denom_tmpl == 0.0) break;
              it_denom_tmpl++;
            }
          }
        }
      }

    }

	}
	if(new_tree) count_trees++;

}

void
coal_LA::init_bootstrap(){

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

}

void 
coal_LA::Dump(std::string filename, std::vector<std::string>& unique_groups){

	//populate num_boot and num_denom (vector of size num_bootstrap)
	//using num, denom, blocks

	init_bootstrap();

	std::vector<std::vector<std::vector<std::vector<double>>>>::iterator it1_num_boot = num_boot.begin(), it1_denom_boot = denom_boot.begin();
	for(it1_blocks = blocks.begin(); it1_blocks != blocks.end(); it1_blocks++){
	
			for(int i = 0; i < num_groups; i++){
				(*it1_num_boot)[i].resize(num_groups);
				for(int j = 0; j < num_groups; j++){
					(*it1_num_boot)[j].resize(num_epochs);
					std::fill((*it1_num_boot)[i][j].begin(), (*it1_num_boot)[i][j].end(), 0.0);
				}
			}
			(*it1_denom_boot).resize(num_groups);
			for(int i = 0; i < num_groups; i++){
				(*it1_denom_boot)[i].resize(num_groups);
				for(int j = 0; j < num_groups; j++){
					(*it1_denom_boot)[j].resize(num_epochs);
					std::fill((*it1_denom_boot)[i][j].begin(), (*it1_denom_boot)[i][j].end(), 0.0);
				}
			}

		it1_num   = num.begin();
		it1_denom = denom.begin();
		assert(num.size() == (*it1_blocks).size());
		for(it2_blocks = (*it1_blocks).begin(); it2_blocks != (*it1_blocks).end(); it2_blocks++){

			if(*it2_blocks > 0){

				assert((*it1_num).size() == (*it1_num_boot).size());

				for(int i = 0; i < num_groups; i++){
					for(int j = 0; j < num_groups; j++){
					  for(int e = 0; e < num_epochs; e++){
						  (*it1_num_boot)[i][j][e]   += *it2_blocks * (*it1_num)[i][j][e];
							(*it1_denom_boot)[i][j][e] += *it2_blocks * (*it1_denom)[i][j][e];
						}
					}
				}

			}

			it1_num++;
			it1_denom++;
		}

		it1_num_boot++;
		it1_denom_boot++;
	}

	std::ofstream os(filename);
	for(int i = 0; i < unique_groups.size(); i++){
		os << unique_groups[i] << " ";
	}
	os << "\n";
	for(it_epochs = epochs.begin(); it_epochs != epochs.end(); it_epochs++){
		os << *it_epochs << " ";
	}
	os << "\n";
	it1_num_boot = num_boot.begin();
	it1_denom_boot = denom_boot.begin();
	for(; it1_num_boot != num_boot.end();){

    for(int i = 0; i < num_groups; i++){
      for(int j = 0; j < num_groups; j++){

				if(i > j){
					os << i << " " << j << " ";
					for(int e = 0; e < num_epochs; e++){
						os << (*it1_num_boot)[i][j][e]/(*it1_denom_boot)[i][j][e] << " ";
					}
					os << "\n";
				}else{
					os << i << " " << j << " ";
					for(int e = 0; e < num_epochs; e++){
						os << (*it1_num_boot)[j][i][e]/(*it1_denom_boot)[j][i][e] << " ";
					}
					os << "\n";
				}
			}
		}

		it1_num_boot++;
		it1_denom_boot++;
	}
	os.close();

}
