#ifndef ADNA_EM_HPP
#define ADNA_EM_HPP

#include "gzstream.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"

#include <iostream>
#include <string>
#include <random>
#include <iterator>


class aDNA_EM{

	private:

		std::vector<double> coal_rates;
		std::vector<double> epochs;
		int num_epochs;

		std::vector<double> A_ep, B_ep, C_ep;

		double logsumexp(double loga, double logb);
		double logminusexp(double loga, double logb);

		void get_tint(double age_begin, double age_end, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin, int& i_end);
		void get_ABC(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B, std::vector<double>& C);
		void get_ABC_lazy(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B, std::vector<double>& C);

	public:

		aDNA_EM(std::vector<double>& epochs, std::vector<double>& coal): epochs(epochs), coal_rates(coal){
			num_epochs = epochs.size();
			A_ep.resize(num_epochs);
			B_ep.resize(num_epochs);
			C_ep.resize(num_epochs);
			std::vector<int> ep_index(num_epochs);
			int i = 0;
			for(std::vector<int>::iterator it_ep_index = ep_index.begin(); it_ep_index != ep_index.end(); it_ep_index++){
				*it_ep_index = i;
				i++;  
			}
			get_ABC(epochs, ep_index, A_ep, B_ep, C_ep);
		}

		double EM_shared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom); 

};

class aDNA_EM_simplified{

	private:

		std::vector<double> coal_rates;
		std::vector<double> epochs;
		int num_epochs;

		std::vector<double> A_ep, B_ep, C_ep;

		double logsumexp(double loga, double logb);
		double logminusexp(double loga, double logb);

		void get_tint(double age, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin);
		void get_ABC(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B);
		void get_ABC_lazy(std::vector<double>& t_int, int i_begin, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B);

	public:

		aDNA_EM_simplified(std::vector<double>& epochs, std::vector<double>& coal): epochs(epochs), coal_rates(coal){
			num_epochs = epochs.size();
			A_ep.resize(num_epochs);
			B_ep.resize(num_epochs);
			std::vector<int> ep_index(num_epochs);
			int i = 0;
			for(std::vector<int>::iterator it_ep_index = ep_index.begin(); it_ep_index != ep_index.end(); it_ep_index++){
				*it_ep_index = i;
				i++;  
			}
			get_ABC(epochs, ep_index, A_ep, B_ep);
		}

		double EM_shared(double age, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared(double age, std::vector<double>& num, std::vector<double>& denom); 

		double EM_shared_exact(double age, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared_exact(double age, std::vector<double>& num, std::vector<double>& denom); 

};

#endif //ADNA_EM_HPP

