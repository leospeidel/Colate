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

class aDNA_EM2{

	private:

		std::vector<double> coal_rates;
		std::vector<double> epochs;
		int num_epochs;

		std::vector<double> A_ep, B_ep;

		double logsumexp(double loga, double logb);
		double logminusexp(double loga, double logb);

		void get_tint(double age_begin, double age_end, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin, int& i_end);
		void get_AB(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B);

	public:

		aDNA_EM2(std::vector<double>& epochs, std::vector<double>& coal): epochs(epochs), coal_rates(coal){
			num_epochs = epochs.size();
			A_ep.resize(num_epochs);
			B_ep.resize(num_epochs);
			std::vector<int> ep_index(num_epochs);
			int i = 0;
			for(std::vector<int>::iterator it_ep_index = ep_index.begin(); it_ep_index != ep_index.end(); it_ep_index++){
				*it_ep_index = i;
				i++;  
			}
			get_AB(epochs, ep_index, A_ep, B_ep);
		}

		double EM_shared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom); 

};

class aDNA_EM_tree{

	private:

		std::vector<double> coal_rates;
		std::vector<double> epochs;
		std::vector<double> t_int;
		std::vector<double> cumsum_coal_rate;
		std::vector<int> ep_index;
		int num_epochs, num_age_bins;
		double C;


		double logsumexp(double loga, double logb);
		double logminusexp(double loga, double logb);

	public:

		aDNA_EM_tree(double C, std::vector<double>& epochs, std::vector<double>& coal): epochs(epochs), coal_rates(coal){
			num_epochs = epochs.size();
			
			num_age_bins = ((int) (log(1e8) * C)) + 1;
			t_int.resize(num_age_bins);
			ep_index.resize(num_age_bins);
			cumsum_coal_rate.resize(num_age_bins);

			std::vector<double>::iterator it_tint = t_int.begin();
			std::vector<int>::iterator it_ep_index = ep_index.begin();
			int bin = 0;
			int e   = 0;
			*it_tint = 0;
			*it_ep_index = 0;
			it_tint++;
			for(bin = 1; bin < num_age_bins; bin++){
				*it_tint =  exp((bin-1.0)/C)/10.0;
				if(e < num_epochs){
					if(e != num_epochs - 1){
				    if(*it_tint > epochs[e+1]) e++;
					}
				}
				*it_ep_index = e;
				it_ep_index++;
				it_tint++;
			}

		}

		void UpdateTree(std::vector<float>& num_lins);
		double EM_shared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom); 

};


class aDNA_EM_tree_fast{

	private:

		std::vector<double> inv_coal_rate_tmpl, sum_coal_rate_tmpl;
		std::vector<double> epochs;
		std::vector<double> t_int;
		std::vector<double> f, tf, tf_prec;
		std::vector<double> cumsum_coal_rate, sum_coal_rate_tint, inv_coal_rate_tint;
		std::vector<int> ep_index;
		std::vector<double> factor;
		int num_epochs, num_age_bins;
		double C;

		//iterators
		std::vector<double>::iterator it_coal;
		std::vector<double>::iterator it_sum_coal;
		std::vector<double>::iterator it_inv_coal;
		std::vector<double>::iterator it_cumsum_coal;
		std::vector<int>::iterator it_ep_index;
		std::vector<double>::iterator it_tint;
		std::vector<double>::iterator it_tint_next;
		std::vector<float>::iterator it_num_lins;
		std::vector<double>::iterator it_f, it_tf, it_tf_prec;

		double logsumexp(double loga, double logb);
		double logminusexp(double loga, double logb);

	public:

		aDNA_EM_tree_fast(double C, std::vector<double>& epochs): epochs(epochs){
			num_epochs = epochs.size();

			num_age_bins = ((int) (log(1e8) * C)) + 1;
			t_int.resize(num_age_bins);
			ep_index.resize(num_age_bins);
			cumsum_coal_rate.resize(num_age_bins);
			sum_coal_rate_tmpl.resize(num_age_bins);
			inv_coal_rate_tmpl.resize(num_age_bins);
			sum_coal_rate_tint.resize(num_age_bins);
			inv_coal_rate_tint.resize(num_age_bins);
			f.resize(num_age_bins);
			tf.resize(num_age_bins);
			tf_prec.resize(num_age_bins);
			factor.resize(num_epochs);

			std::vector<double>::iterator it_tint = t_int.begin();
			std::vector<int>::iterator it_ep_index = ep_index.begin();
			int bin = 0;
			int e   = 0;
			*it_tint = 0;
			*it_ep_index = 0;
			it_tint++;
			for(bin = 1; bin < num_age_bins; bin++){
				*it_tint =  exp((bin-1.0)/C)/10.0;
				if(e < num_epochs){
					if(e != num_epochs - 1){
						if(*it_tint > epochs[e+1]) e++;
					}
				}
				*it_ep_index = e;
				it_ep_index++;
				it_tint++;
			}

		}

		void UpdateCoal(std::vector<double>& icoal);
		void UpdateTree(std::vector<float>& num_lins);
		double EM_shared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom);
		double EM_notshared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom); 

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

		double EM_shared_exact(double age, std::vector<double>& num, std::vector<double>& denom, double coal_rate);
		double EM_notshared_exact(double age, std::vector<double>& num, std::vector<double>& denom, double coal_rate); 

};

#endif //ADNA_EM_HPP

