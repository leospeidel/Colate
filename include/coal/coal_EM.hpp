#ifndef COAL_EM_HPP
#define COAL_EM_HPP

#include "gzstream.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"

#include <iostream>
#include <string>
#include <random>
#include <iterator>

class coal_EM{

	private:

    std::mt19937 rng;
    int seed = 2;

    double log_0 = std::log(0.0);

		std::vector<double> coal_rates;
		std::vector<double> epochs;
		std::vector<int> ep_index_tmpl;
		int num_epochs;

		std::vector<double> A_ep, B_ep;

		double logsumexp(double loga, double logb);
		double logminusexp(double loga, double logb);

		void get_tint(double age_begin, double age_end, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin, int& i_end);
		void get_AB(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B);

	public:

		coal_EM(std::vector<double>& epochs, std::vector<double>& coal): epochs(epochs), coal_rates(coal){  
      rng.seed(seed);
			num_epochs = epochs.size();
			A_ep.resize(num_epochs);
			B_ep.resize(num_epochs);
			ep_index_tmpl.resize(num_epochs);
			int i = 0;
			for(std::vector<int>::iterator it_ep_index = ep_index_tmpl.begin(); it_ep_index != ep_index_tmpl.end(); it_ep_index++){
				*it_ep_index = i;
				i++;  
			}
			get_AB(epochs, ep_index_tmpl, A_ep, B_ep);
		}

		void UpdateCoal(std::vector<double>& coal){
      coal_rates = coal;
			get_AB(epochs, ep_index_tmpl, A_ep, B_ep);
		}

		double EM_shared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom); 

		double EM_shared_sampled(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared_sampled(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom); 

};

#endif //COAL_EM_HPP

