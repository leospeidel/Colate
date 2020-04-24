#include "catch.hpp"

#include "data.hpp"
#include "anc.hpp"
#include "coal_EM.hpp"

double
logsumexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isnan(loga)) std::cerr << loga << std::endl;
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return logb;
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isnan(logb)) std::cerr << logb << std::endl;
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	if(loga > logb){
		return(loga + log(1.0 + exp(logb - loga)));
	}else{
		return(logb + log(1.0 + exp(loga - logb)));
	}

}

double
logminusexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isnan(loga)) std::cerr << loga << std::endl;
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return log(0.0); //assuming small value
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isnan(logb)) std::cerr << logb << std::endl;
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	//assert(loga > logb);
	if(loga < logb){
		//std::cerr << loga << " " << logb << std::endl;
		return(log(0));
	}
	return(loga + log1p(-exp(logb - loga)));
	//return(loga + log(-expm1(logb-loga)));
	//return(loga + log(1.0 - exp(logb - loga)));

}


TEST_CASE("test EM expectation step"){

	double p = 1e-1;

	//decide on epochs
	float years_per_gen = 28.0;
	int num_epochs = 20;
	num_epochs++;
	std::vector<double> epochs(num_epochs);
	epochs[0] = 0.0;
	epochs[1] = 1e3/years_per_gen;
	float log_10 = std::log(10);
	for(int e = 2; e < num_epochs-1; e++){
		epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
	}
	epochs[num_epochs-1] = 1e8/years_per_gen;

	double coal_rate = 1e-7;

	for(int f = 1; f <= 7; f++){

		double initial_coal_rate = coal_rate * exp(log(10) * (f-1));
    std::cerr << initial_coal_rate << std::endl;

		Data data(1, 1);
		std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);
		for(int i = 0; i < data.N; i++){
			coal_rates[i].resize(num_epochs);
			coal_rates_num[i].resize(num_epochs);
			coal_rates_denom[i].resize(num_epochs);
			std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
		}

		coal_EM EM(epochs, coal_rates[0]);
		coal_EM_simplified EM2(epochs, coal_rates[0]);

		int i = 0;
		std::vector<double> num(num_epochs), denom(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1
		std::vector<double> num2(num_epochs), denom2(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1

		double C = 5;
		int num_age_bins = ((int) (log(1e8) * C));
		std::vector<int> age_shared_count(num_age_bins, 0), age_notshared_count(num_age_bins, 0);
		std::vector<double> age_bin(num_age_bins, 0.0);
		int bin = 0;
		for(std::vector<double>::iterator it_age_bin = age_bin.begin(); it_age_bin != age_bin.end(); it_age_bin++){
			*it_age_bin = exp(bin/C)/10.0;
			bin++;
		}

		//shared
		double lambda = initial_coal_rate;
		for(int bin1 = 0; bin1 < num_age_bins; bin1++){
			int bin2 = bin1;
			std::fill(num.begin(), num.end(), 0.0);
			std::fill(denom.begin(), denom.end(), 0.0);
			double logl = EM.EM_shared(age_bin[bin1], age_bin[bin2], num, denom);

			std::fill(num2.begin(), num2.end(), 0.0);
			std::fill(denom2.begin(), denom2.end(), 0.0);
			//double logl2 = EM2.EM_shared_exact(age_bin[bin1], num2, denom2, initial_coal_rate);
			double logl2 = EM2.EM_shared(age_bin[bin1], num2, denom2);

			REQUIRE(std::fabs(logl - logl2) < 1e-3);
			//std::cerr << "logl " << logl << " " << logl2 << std::endl;
			for(int e = 0; e < num_epochs-1; e++){
				if(1){	
					//std::cerr << e << " " << num[e] << " " << denom[e] << std::endl;
					//std::cerr << e << " " << num2[e] << " " << denom2[e] << std::endl;
					//std::cerr << e << " " << num[e] - num2[e] << " " << denom[e] - denom2[e] << std::endl << std::endl;
          if(std::fabs(denom[e] - denom2[e]) > p) exit(1);
					if(num2[e] > 0){
						REQUIRE( (std::fabs(num[e] - num2[e]) <= p || std::fabs(num[e] - num2[e])/num2[e] <= p) );
					}else{
						REQUIRE(std::fabs(num[e]) < p);
					}
					if(denom2[e] > 0){
						REQUIRE( (std::fabs(denom[e] - denom2[e]) <= p || std::fabs(denom[e] - denom2[e])/denom2[e] <= p) );
					}else{
						REQUIRE(std::fabs(denom[e]) < p);
					}
				}
			}

		}

		for(int bin1 = 0; bin1 < num_age_bins; bin1++){
			int bin2 = bin1;
			std::fill(num.begin(), num.end(), 0.0);
			std::fill(denom.begin(), denom.end(), 0.0);
			double logl = EM.EM_notshared(age_bin[bin1], age_bin[bin2], num, denom);

			std::fill(num2.begin(), num2.end(), 0.0);
			std::fill(denom2.begin(), denom2.end(), 0.0);
			//double logl2 = EM2.EM_notshared_exact(age_bin[bin1], num2, denom2, initial_coal_rate);
			double logl2 = EM2.EM_notshared(age_bin[bin1], num2, denom2);

			REQUIRE(std::fabs(logl - logl2) < p);

			for(int e = 0; e < num_epochs-1; e++){
				if(1){	
					//std::cerr << e << " " << num[e] << " " << denom[e] << std::endl;
					//std::cerr << e << " " << num2[e] << " " << denom2[e] << std::endl;
					//std::cerr << e << " " << num[e] - num2[e] << " " << denom[e] - denom2[e] << std::endl << std::endl;
					if(num2[e] > 0){
						REQUIRE( (std::fabs(num[e] - num2[e]) <= p || std::fabs(num[e] - num2[e])/num2[e] <= p) );
					}else{
						REQUIRE(std::fabs(num[e]) < p);
					}
					if(denom2[e] > 0){
						REQUIRE( (std::fabs(denom[e] - denom2[e]) <= p || std::fabs(denom[e] - denom2[e])/denom2[e] <= p) );
					}else{
						REQUIRE(std::fabs(denom[e]) < p);
					}
				}
			}

		}

		for(int bin1 = 0; bin1 < num_age_bins; bin1++){
			for(int bin2 = bin1; bin2 < num_age_bins; bin2++){
				std::fill(num.begin(), num.end(), 0.0);
				std::fill(denom.begin(), denom.end(), 0.0);
				EM.EM_shared(age_bin[bin1], age_bin[bin2], num, denom);
				for(int e = 0; e < num_epochs; e++){
					REQUIRE(!std::isnan(num[e]));
					REQUIRE(!std::isnan(denom[e]));
					REQUIRE(num[e] >= 0.0);
					REQUIRE(denom[e] >= 0.0);
				}
				std::fill(num.begin(), num.end(), 0.0);
				std::fill(denom.begin(), denom.end(), 0.0);
				EM.EM_notshared(age_bin[bin1], age_bin[bin2], num, denom);
				for(int e = 0; e < num_epochs; e++){
					REQUIRE(!std::isnan(num[e]));
					REQUIRE(!std::isnan(denom[e]));
					REQUIRE(num[e] >= 0.0);
					REQUIRE(denom[e] >= 0.0);
				}
			}
		}

	}

}


TEST_CASE("test EM_tree expectation step"){

	double p = 1e-1;

	//decide on epochs
	float years_per_gen = 28.0;
	int num_epochs = 20;
	num_epochs++;
	std::vector<double> epochs(num_epochs);
	epochs[0] = 0.0;
	epochs[1] = 1e3/years_per_gen;
	float log_10 = std::log(10);
	for(int e = 2; e < num_epochs-1; e++){
		epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
	}
	epochs[num_epochs-1] = 1e8/years_per_gen;

	double coal_rate = 1e-7;

	for(int f = 1; f <= 6; f++){

		double initial_coal_rate = coal_rate * exp(log(10) * (f-1));
		std::cerr << initial_coal_rate << std::endl;

		Data data(1, 1);
		std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);
		for(int i = 0; i < data.N; i++){
			coal_rates[i].resize(num_epochs);
			coal_rates_num[i].resize(num_epochs);
			coal_rates_denom[i].resize(num_epochs);
			std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
		}
	
		int i = 0;
		std::vector<double> num(num_epochs), denom(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1
		std::vector<double> num2(num_epochs), denom2(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1
		std::vector<double> num3(num_epochs), denom3(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1

		double C = 5;
		int num_age_bins = ((int) (log(1e8) * C)) + 1;
		std::vector<int> age_shared_count(num_age_bins, 0), age_notshared_count(num_age_bins, 0);
		std::vector<double> age_bin(num_age_bins, 0.0);
		int bin = 1;
		for(std::vector<double>::iterator it_age_bin = std::next(age_bin.begin(),1); it_age_bin != age_bin.end(); it_age_bin++){
			*it_age_bin = exp((bin-1.0)/C)/10.0;
			bin++;
		}

		std::vector<float> DAF(num_age_bins, 2), DAF2(num_age_bins, 1), num_lins(num_age_bins, 2);

		coal_EM_tree_fast EM(C, 2, epochs);
		EM.UpdateCoal(coal_rates[0]);
		coal_EM_tree EM3(C, epochs, coal_rates[0]);
		EM.UpdateTree(num_lins);
		for(int e = 0; e < num_epochs; e++){
      coal_rates[0][e] *= 2;
		}
		//coal_EM_simplified EM2(epochs, coal_rates[0]);
		coal_EM EM2(epochs, coal_rates[0]);
    EM.UpdateTree(num_lins);
		EM3.UpdateTree(num_lins);

		//shared
		double lambda = initial_coal_rate;
		for(int bin1 = 1; bin1 < num_age_bins; bin1++){
			int bin2 = bin1;

			for(int i = 0; i < bin1; i++){
        DAF[i] = 2;
				DAF2[i] = 1;
			}
			for(int i = bin1; i < num_age_bins; i++){
        DAF[i] = 0;
				DAF2[i] = 0;
			}

			double logl, logl2, logl3, logl4;

			std::fill(num.begin(), num.end(), 0.0);
			std::fill(denom.begin(), denom.end(), 0.0);
			logl = EM.EM_shared(age_bin[bin1], age_bin[bin2], num_lins, DAF2, num, denom);

			std::fill(num2.begin(), num2.end(), 0.0);
			std::fill(denom2.begin(), denom2.end(), 0.0);
			logl2 = EM2.EM_shared(age_bin[bin1], age_bin[bin2], num2, denom2);

			std::fill(num3.begin(), num3.end(), 0.0);
			std::fill(denom3.begin(), denom3.end(), 0.0);
			logl3 = EM3.EM_shared(age_bin[bin1], age_bin[bin2], num_lins, DAF2, num2, denom2);

			logl4 = EM.Logl_shared(age_bin[bin1], age_bin[bin2], num_lins, DAF2);

			REQUIRE(std::fabs(logl - logl2) < 1e-3);
			REQUIRE(std::fabs(logl - logl3) < 1e-3);
			REQUIRE(std::fabs(logl - logl4) < 1e-3);

			/*
			std::cerr << bin1 << " " << logl << " " << logl3 << std::endl;
			for(int e = 0; e < num_epochs; e++){
		  	std::cerr << e << " " << num[e] << " " << num2[e] << " " << denom[e] << " " << denom2[e] << std::endl;
        REQUIRE(std::fabs(num[e] - num2[e]) < 1e-3);
				REQUIRE(std::fabs(denom[e] - denom2[e]) < 1e-3);
			}
			*/

			if(bin1 < num_age_bins-1){
			std::fill(num.begin(), num.end(), 0.0);
			std::fill(denom.begin(), denom.end(), 0.0);
			logl = EM.EM_notshared(age_bin[bin1], age_bin[bin2], num_lins, DAF2, num, denom);

			std::fill(num2.begin(), num2.end(), 0.0);
			std::fill(denom2.begin(), denom2.end(), 0.0);
			logl2 = EM2.EM_notshared(age_bin[bin1], age_bin[bin2], num2, denom2);

			std::fill(num3.begin(), num3.end(), 0.0);
			std::fill(denom3.begin(), denom3.end(), 0.0);
			logl3 = EM3.EM_notshared(age_bin[bin1], age_bin[bin2], num_lins, DAF2, num2, denom2);

			logl4 = EM.Logl_notshared(age_bin[bin1], age_bin[bin2], num_lins, DAF2);

			if(!std::isinf(logl)){
				//std::cerr << bin1 << " " << logl << " " << logl3 << std::endl;
				REQUIRE(std::fabs(logl - logl3) < 1e-3);
				REQUIRE(std::fabs(logl - logl4) < 1e-3);
				for(int e = 0; e < num_epochs-1; e++){
					//std::cerr << e << " " << num[e] << " " << num2[e] << " " << denom[e] << " " << denom2[e] << std::endl;
					REQUIRE(std::fabs(num[e] - num2[e]) < 1e-2);
					REQUIRE(std::fabs(denom[e] - denom2[e]) < 1e-2);
				}
			}

			}

		}
	
	}

}


TEST_CASE("test simplified EM expectation step"){

	//decide on epochs
	float years_per_gen = 28.0;
	int num_epochs = 30;
	num_epochs++;
	std::vector<double> epochs(num_epochs);
	epochs[0] = 0.0;
	epochs[1] = 1e3/years_per_gen;
	float log_10 = std::log(10);
	for(int e = 2; e < num_epochs-1; e++){
		epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
	}
	epochs[num_epochs-1] = 1e8/years_per_gen;

	double initial_coal_rate = 1.0/5000.0;
	Data data(1, 1);
	std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);
	for(int i = 0; i < data.N; i++){
		coal_rates[i].resize(num_epochs);
		coal_rates_num[i].resize(num_epochs);
		coal_rates_denom[i].resize(num_epochs);
		std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
	}

	coal_EM_simplified EM(epochs, coal_rates[0]);

	int i = 0;
	std::vector<double> num(num_epochs), denom(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1

	double C = 5;
	int num_age_bins = ((int) (log(1e8) * C));
	std::cerr << num_age_bins << std::endl;
	std::vector<int> age_shared_count(num_age_bins, 0), age_notshared_count(num_age_bins, 0);
	std::vector<double> age_bin(num_age_bins, 0.0);
	int bin = 0;
	for(std::vector<double>::iterator it_age_bin = age_bin.begin(); it_age_bin != age_bin.end(); it_age_bin++){
		*it_age_bin = exp(bin/C)/10.0;
		bin++;
	}

	double lambda = initial_coal_rate;
	for(int bin = 0; bin < num_age_bins; bin++){

		std::fill(num.begin(), num.end(), 0.0);
		std::fill(denom.begin(), denom.end(), 0.0);
		double logl = EM.EM_shared(age_bin[bin], num, denom);

		for(int e = 0; e < num_epochs-1; e++){

			double ep1 = epochs[e], ep2 = epochs[e+1];

			double A = std::min(ep1, age_bin[bin]), B = std::min(age_bin[bin], ep2);
			double norm     = logminusexp(0.0, -lambda*age_bin[bin]);
			double th_num   = logminusexp(-lambda*A, -lambda*B);

			/*
				 double th_denom = logsumexp(log(A)-lambda*A, th_num - log(lambda));
				 th_denom = logminusexp(th_denom, log(B) - lambda*B);
				 th_denom = logminusexp(th_denom, log(ep1) + th_num);
				 */

			double th_denom = log(A*exp(-lambda*A) + exp(th_num)/lambda - B*exp(-lambda*B) - ep1*exp(th_num));

			double rest = logminusexp(-lambda*B, -lambda*age_bin[bin]);
			th_denom = logsumexp(th_denom, log(ep2 - ep1) + rest);

			th_num   -= norm;
			th_denom -= norm;

			th_num = exp(th_num);
			th_denom = exp(th_denom);

			if(1){	
				REQUIRE(std::fabs(logl - norm) < 1e-3);
				if(th_num > 0){
					REQUIRE(std::fabs(num[e] - th_num) <= 1e-5);
				}else{
					REQUIRE(num[e] == 0.0);
				}
				if(th_denom > 0){
					REQUIRE(std::fabs(denom[e] - th_denom) <= 1e-5);
				}else{
					REQUIRE(denom[e] == 0.0);
				}
			}

		}
	}

	for(int bin = 0; bin < num_age_bins; bin++){

		std::fill(num.begin(), num.end(), 0.0);
		std::fill(denom.begin(), denom.end(), 0.0);
		double logl = EM.EM_notshared(age_bin[bin], num, denom);

		for(int e = 0; e < num_epochs-1; e++){

			double ep1 = epochs[e], ep2 = epochs[e+1];

			double A = std::max(ep1, age_bin[bin]), B = std::max(age_bin[bin], ep2);
			double norm     = -lambda*age_bin[bin];
			double th_num   = logminusexp(-lambda*A, -lambda*B);

			/*
				 double th_denom = logsumexp(log(A)-lambda*A, th_num - log(lambda));
				 th_denom = logminusexp(th_denom, log(B) - lambda*B);
				 th_denom = logminusexp(th_denom, log(ep1) + th_num);
				 */
			double th_denom = log(A*exp(-lambda*A) + exp(th_num)/lambda - B*exp(-lambda*B) - ep1*exp(th_num));

			double rest = -lambda*B;
			th_denom = logsumexp(th_denom, log(ep2 - ep1) + rest);

			th_num   -= norm;
			th_denom -= norm;

			th_num = exp(th_num);
			th_denom = exp(th_denom);

			if(1){
				REQUIRE(std::fabs(logl - norm) < 1e-3);
				if(th_num > 0){
					REQUIRE(std::fabs(num[e] - th_num) <= 1e-5);
				}else{
					REQUIRE(num[e] == 0.0);
				}
				if(th_denom > 0){
					REQUIRE(std::fabs(denom[e] - th_denom) <= 1e-5);
				}else{
					REQUIRE(denom[e] == 0.0);
				}
			}

		}
	}


}

