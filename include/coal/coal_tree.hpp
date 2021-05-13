#ifndef COAL_TREE_EM_HPP
#define COAL_TREE_EM_HPP

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <random>
#include <ctime>

#include "anc.hpp"
#include "mutations.hpp"
#include "cxxopts.hpp"

//input tree, compute MLE of coalescence rates for tree

class coal_tree {

  private:

		std::mt19937 rng;

    int N, N_total, num_trees, block_size, num_blocks, num_bootstrap, current_block, count_trees, num_epochs;
		std::vector<float> coords;
		std::vector<float>::iterator it_coords, it_coords_next;
		std::vector<int> num_lins, sorted_indices;
		std::vector<int>::iterator it_num_lins, it_sorted_indices, it_sorted_indices_prev;
    std::vector<double> epochs;
		std::vector<double>::iterator it_epochs;

		std::vector<std::vector<int>> blocks;
		std::vector<std::vector<int>>::iterator it1_blocks;
		std::vector<int>::iterator it2_blocks;

		std::vector<std::vector<double>> num, denom, num_boot, denom_boot;
		std::vector<std::vector<double>>::iterator it1_num, it1_denom;
		std::vector<double>::iterator it2_num, it2_denom;

	public:

		coal_tree(std::vector<double> epochs, int num_bootstrap, int block_size);
		coal_tree(std::vector<double> epochs, int num_bootstrap, int block_size, AncMutIterators& ancmut);

    void update_ancmut(AncMutIterators& ancmut);
		void populate(Tree& tree, double num_bases_tree_persists);
    void init_bootstrap();
		void Dump(std::string filename);

};

class coal_LA {

	private:

		std::mt19937 rng;

		int N, N_total, num_trees, block_size, num_blocks, num_bootstrap, num_groups, current_block, count_trees, num_epochs;
		std::vector<Leaves> desc;
		std::vector<float> coords;
		std::vector<int> sorted_indices;
		std::vector<int>::iterator it_sorted_indices;
		std::vector<double> epochs;
		std::vector<double>::iterator it_epochs;

		std::vector<std::vector<int>> blocks;
		std::vector<std::vector<int>>::iterator it1_blocks;
		std::vector<int>::iterator it2_blocks;

		std::vector<std::vector<std::vector<std::vector<double>>>> num, denom, num_boot, denom_boot;
		std::vector<std::vector<std::vector<std::vector<double>>>>::iterator it1_num, it1_denom;

	public:

		coal_LA(std::vector<double> epochs, int num_bootstrap, int block_size, int num_groups);
		coal_LA(std::vector<double> epochs, int num_bootstrap, int block_size, int num_groups, AncMutIterators& ancmut);

		void update_ancmut(AncMutIterators& ancmut);
		void populate(Tree& tree, double num_bases_tree_persists, std::vector<int>& group, bool new_tree);
		void init_bootstrap();
		void Dump(std::string filename, std::vector<std::string>& unique_groups);

};

#endif //COAL_TREE_HPP

