#include "coal_EM.hpp"


////////////////////////////////////

double
coal_EM::logsumexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return logb;
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	if(loga > logb){
		long double res = loga + log1p(exp(logb - loga));
		return(res);
	}else{
		long double res = logb + log1p(exp(loga - logb));
		return(res);
	}

}

double
coal_EM::logminusexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return log(0.0); //assuming small value
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	if(loga < logb){
		return(log(0));
	}
	return(loga + log1p(-exp(logb - loga)));
	//return(loga + log(-expm1(logb-loga)));
	//return(loga + log(1.0 - exp(logb - loga)));

}

void
coal_EM::get_tint(double age_begin, double age_end, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin, int& i_end){

	int e = 0;
	int i = 0;
	for(e = 0; e < num_epochs; e++){
		if(age_begin < epochs[e] && i_begin == -1){
			t_int[i] = age_begin;
			ep_index[i] = e-1;
			i_begin = i;
			i++;
		}
		if(age_end   < epochs[e] && i_end == -1){
			t_int[i] = age_end;
			ep_index[i] = e-1;
			i_end = i;
			i++;
		}
		t_int[i] = epochs[e];
		ep_index[i] = e;
		i++;
	}
	if(i_begin == -1){
		t_int[i] = age_begin;
		ep_index[i] = num_epochs-1;
		i_begin = i;
		i++;
	}
	if(i_end == -1){
		t_int[i] = age_end;
		ep_index[i] = num_epochs-1;
		i_end = i;
		i++;
	}

}

void 
coal_EM::get_AB(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B){

	std::vector<double> cumsum_coal_rate(t_int.size(), 0.0);
	for(int i = 1; i < t_int.size(); i++){
		cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * (t_int[i] - t_int[i-1]);
	}

	double t_begin, t_end, coal_rate_e, inv_coal_rate_e;
	int i = 0;
	for(i = 0; i < t_int.size()-1; i++){

		t_begin     = t_int[i];
		t_end       = t_int[i+1];
		coal_rate_e = coal_rates[ep_index[i]];
		inv_coal_rate_e = 1.0/coal_rates[ep_index[i]];

		assert(t_end >= t_begin);
		assert(ep_index[i] == i);

		if(coal_rate_e > 0 && t_end != 0 && t_end - t_begin > 0){

			A[i] = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
			B[i] = (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]);
			B[i] = log(B[i]) - cumsum_coal_rate[i];

		}else{

			A[i] = log(0.0); 
			B[i] = log(0.0);

		}

		assert(!std::isnan(A[i]));
		assert(!std::isnan(B[i]));

	}
	i = t_int.size()-1;
	t_begin = t_int[i];
	coal_rate_e = coal_rates[ep_index[i]];
	if(coal_rate_e > 0){

		A[i] = -cumsum_coal_rate[i];
		B[i] = log(t_begin + 1.0/coal_rate_e) - cumsum_coal_rate[i];

	}else{

		A[i] = log(0.0); 
		B[i] = log(0.0);

	}
	assert(!std::isnan(A[i]));
	assert(!std::isnan(B[i]));

}

double
coal_EM::EM_shared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom){

	//calculate A, B, then f, and tf
	std::fill(num.begin(), num.end(), 0);
	std::fill(denom.begin(), denom.end(), 0);

	bool times_identical = false;
	if(age_begin == age_end) times_identical = true;

	int num_epochs = epochs.size();

	assert(age_begin <= age_end);
	bool age_identical = false;
	if(age_begin == age_end) age_identical = true;

	//make vector of time intervals
	std::vector<double> t_int(num_epochs + 2, 0.0);
	std::vector<int> ep_index(num_epochs + 2, 0);
	int i_begin = -1, i_end = -1; 
	get_tint(age_begin, age_end, t_int, ep_index, i_begin, i_end);
	assert(ep_index[i_end] == ep_index[i_end - 1]);

	std::vector<double> cumsum_coal_rate(t_int.size(), 0.0);
	for(int i = 1; i < t_int.size(); i++){
		cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * (t_int[i] - t_int[i-1]);
	}

	double t_begin, t_end, coal_rate_e, inv_coal_rate_e;
	int current_e = -1;
	double normalising_constant = 1.0;

	int i = 0;
	for(int e = 0; e < ep_index[i_end] + 1; e++){

		coal_rate_e = coal_rates[e];
		inv_coal_rate_e = 1.0/coal_rates[e];
		double num_e, denom_e;
		if(e < ep_index[i_begin]){

			assert(i == e);
			num[e]   = A_ep[i];
			denom[e] = B_ep[i];
			i++;

		}else if(e == ep_index[i_begin]){
			assert(i+1 == i_begin);
			t_begin  = t_int[i];
			t_end    = t_int[i+1];
			num_e   = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
			denom_e = log( (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e)*exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) ) - cumsum_coal_rate[i];
			i++;
		}

		if(e >= ep_index[i_begin] && !times_identical){
			t_begin = t_int[i];
			t_end   = t_int[i+1];

			num[e] = log((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i])) - cumsum_coal_rate[i] - log(age_end - age_begin);

			double x_begin = t_begin/age_end, x_end = t_end/age_end;
			double term1 = (x_begin*(age_end-t_begin)/inv_coal_rate_e + 1.0 - 2.0*(x_begin + inv_coal_rate_e/age_end));
			double term2 = (-x_end*(age_end -t_end)/inv_coal_rate_e   - 1.0 + 2.0*(x_end   + inv_coal_rate_e/age_end));
			double tmp = term1;
			tmp       += exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) * term2;

			if(tmp < 0.0){
				denom[e] = log(0.0);
			}else{
				tmp        = log(tmp);
				tmp       += log(age_end) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
				denom[e]   = tmp - log(age_end-age_begin);
			}

			i++;
			while(ep_index[i] == e){
				i++;
				if(i == ep_index.size()) break;
			}
		}

		if(e == ep_index[i_begin]){
			if(times_identical){
				num[e]   = num_e;
				denom[e] = denom_e;
			}else{
				num[e]   = logsumexp(num[e], num_e);
				denom[e] = logsumexp(denom[e], denom_e);
			}
		}

		if(normalising_constant == 1.0){
			normalising_constant = num[e];
		}else{
			normalising_constant = logsumexp(normalising_constant, num[e]);
		}

		if(i == ep_index.size()) break;
	}

	double integ = 1.0;
	int e = 0;
	for(e = 0; e < std::min(num_epochs-1, ep_index[i_end] + 1); e++){
		num[e]   -= normalising_constant;
		denom[e] -= normalising_constant;
		num[e] = exp(num[e]);
		if(integ > 0.0){
			integ -= num[e];
		}else{
			integ  = 0.0;
		}
		denom[e]  = exp(denom[e]);
		denom[e] += -epochs[e] * num[e] + (epochs[e+1]-epochs[e])*integ;
		if(denom[e] < 0.0) denom[e] = 0.0;
	}
	if(ep_index[i_end] == num_epochs - 1){
		e = num_epochs-1;
		num[e]   -= normalising_constant;
		denom[e] -= normalising_constant;
		num[e]    = exp(num[e]);
		denom[e]  = exp(denom[e]);
		denom[e] -= epochs[e]*num[e];
		if(denom[e] < 0.0) denom[e] = 0.0;
	}

	return normalising_constant;

}

double
coal_EM::EM_notshared(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom){

	//calculate A, B, C, then f, and tf
	int num_epochs = epochs.size();

	assert(age_begin <= age_end);
	bool times_identical = false;
	if(age_begin == age_end) times_identical = true;

	//make vector of time intervals
	std::vector<double> t_int(num_epochs + 2, 0.0);
	std::vector<int> ep_index(num_epochs + 2, 0);
	int i_begin = -1, i_end = -1; 
	get_tint(age_begin, age_end, t_int, ep_index, i_begin, i_end);
	assert(ep_index[i_end] == ep_index[i_end - 1]);

	std::vector<double> cumsum_coal_rate(t_int.size(), 0.0);
	for(int i = 1; i < t_int.size(); i++){
		cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * (t_int[i] - t_int[i-1]);
	}

	double t_begin, t_end, coal_rate_e, inv_coal_rate_e;
	int current_e = -1;
	double normalising_constant = 1.0;

	int i = i_begin;
	int e = ep_index[i];
	coal_rate_e = coal_rates[e];
	inv_coal_rate_e = 1.0/coal_rates[e];
	if(times_identical){

		assert(i+1 == i_end);
		i++;
		if(e != num_epochs-1){
			t_begin  = t_int[i];
			t_end    = t_int[i+1];
			num[e]   = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
			denom[e] = log( (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) ) - cumsum_coal_rate[i];
			normalising_constant = num[e];
			i++;
			e++;
			for(; e < num_epochs; e++){
				num[e]   = A_ep[e];
				denom[e] = B_ep[e];
				normalising_constant = logsumexp(normalising_constant, num[e]);
			}
		}else{
			t_begin   = t_int[i];
			num[e]    = -cumsum_coal_rate[i];
			denom[e]  = log(t_begin + inv_coal_rate_e)-cumsum_coal_rate[i];
			normalising_constant = num[e];
			i++;
		}

	}else{

		for(; e < ep_index[i_end]+1; e++){

			coal_rate_e = coal_rates[e];
			inv_coal_rate_e = 1.0/coal_rates[e];
			assert(ep_index[i] == e);
			t_begin = t_int[i];
			t_end   = t_int[i+1];

			num[e] = log( (t_begin - age_begin + inv_coal_rate_e) + (age_begin - t_end - inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i])) - cumsum_coal_rate[i] - log(age_end - age_begin);

			double x_begin = t_begin/age_end, x_end = t_end/age_end, x_age_begin = age_begin/age_end;

			double term1 = (x_begin*(t_begin - age_begin)/inv_coal_rate_e + 2.0*(x_begin + inv_coal_rate_e/age_end) - x_age_begin);
			double term2 = (-x_end*(t_end - age_begin)/inv_coal_rate_e    - 2.0*(x_end   + inv_coal_rate_e/age_end) + x_age_begin);

			double tmp = term1;
			tmp       += exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) * term2;

			if(tmp < 0.0){
				denom[e] = log(0.0);
			}else{
				tmp        = log(tmp);
				tmp       += log(age_end) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
				denom[e]   = tmp - log(age_end - age_begin);
			}
			i++;
			if(ep_index[i_end] == e){

				assert(i == i_end);
				if(e != num_epochs-1){
					t_begin  = t_int[i];
					t_end    = t_int[i+1];
					num[e]   = logsumexp(num[e], logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]));
					denom[e] = logsumexp(denom[e], log( (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i]) ) - cumsum_coal_rate[i]);
					i++;
				}else{
					t_begin   = t_int[i];
					num[e]    = logsumexp(num[e], -cumsum_coal_rate[i]);
					denom[e]  = logsumexp(denom[e], log(t_begin + inv_coal_rate_e)-cumsum_coal_rate[i]);
					i++;
				}

			}

			if(normalising_constant == 1.0){
				normalising_constant = num[e];
			}else{
				normalising_constant = logsumexp(normalising_constant, num[e]);
			}
		}

		for(; e < num_epochs; e++){
			num[e] = A_ep[e];
			denom[e] = B_ep[e];
			normalising_constant = logsumexp(normalising_constant, num[e]);
		}

	}

	double integ = 1.0;
	for(e = 0; e < ep_index[i_begin]; e++){
		num[e]   = 0.0;
		denom[e] = (epochs[e+1]-epochs[e]);
	}
	for(; e < num_epochs-1; e++){
		num[e]   -= normalising_constant;
		denom[e] -= normalising_constant;
		num[e] = exp(num[e]);
		if(integ > 0.0){
			integ -= num[e];
		}else{
			integ  = 0.0;
		}
		denom[e]  = exp(denom[e]);
		denom[e] += -epochs[e] * num[e] + (epochs[e+1]-epochs[e])*integ;
		if(denom[e] < 0.0) denom[e] = 0.0;
	}
	e = num_epochs-1;
	num[e]   -= normalising_constant;
	denom[e] -= normalising_constant;
	num[e]    = exp(num[e]);
	denom[e]  = exp(denom[e]);
	denom[e] -= epochs[e]*num[e];
	if(denom[e] < 0.0) denom[e] = 0.0;
	return normalising_constant;

}


////////////////////////////////////

double
coal_EM_tree::logsumexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return logb;
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	if(loga > logb){
		long double res = loga + log1p(exp(logb - loga));
		return(res);
	}else{
		long double res = logb + log1p(exp(loga - logb));
		return(res);
	}

}

double
coal_EM_tree::logminusexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return log(0.0); //assuming small value
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	if(loga < logb){
		return(log(0));
	}
	return(loga + log1p(-exp(logb - loga)));
	//return(loga + log(-expm1(logb-loga)));
	//return(loga + log(1.0 - exp(logb - loga)));

}

void
coal_EM_tree::UpdateTree(std::vector<float>& num_lins){

	for(int i = 1; i < t_int.size(); i++){
		cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * num_lins[i-1] * (t_int[i] - t_int[i-1]);
	}

}

double
coal_EM_tree::EM_shared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom){

	//prespecify t_int as a time grid, and ep_index according to t_int
	//input num_lin and DAF according to t_int
	//calculate integrals and sum up
	std::fill(num.begin(), num.end(), log(0.0));
	std::fill(denom.begin(), denom.end(), log(0.0));

	assert(age_begin <= age_end);
	bool age_identical = false;
	if(age_begin == age_end) age_identical = true;

	assert(num_lins.size() == num_age_bins);
	assert(DAF.size() == num_age_bins);

	double t_begin, t_end, inv_coal_rate_e;
	int current_e = -1;
	double normalising_constant = 1.0;

	int i = 0, e = 0;
	for(e = 0; e < num_epochs; e++){

		double num_e, denom_e;
		double cumsum_bl = 0.0;

		while(ep_index[i] == e){

			inv_coal_rate_e = 1.0/(num_lins[i] * coal_rates[e]);
			t_begin         = t_int[i];
			t_end           = t_int[i+1];

			if(t_end <= age_begin){
				//constant regime
				num_e   = log(DAF[i]) + logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
				assert(!std::isnan(num_e));
				num[e]  = logsumexp(num[e], num_e);

				denom_e  = log(DAF[i]) + log((t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i])) - cumsum_coal_rate[i];
				assert(!std::isnan(denom_e));
				denom_e  = logminusexp(denom_e, log(t_begin) + num_e) + log(num_lins[i]);
				assert(!std::isnan(denom_e));
				denom[e] = logsumexp(denom[e], logsumexp(denom_e, log(cumsum_bl) + num_e));

				//if(e == 1) std::cerr << num_e << " ";
				//if(e == 0) std::cerr << logsumexp(denom_e, log(cumsum_bl) + num_e) << " ";

			}else if(t_begin >= age_begin && t_end <= age_end && !age_identical){

				if(1){

					//linear regime
					num_e  = log(DAF[i]) + 
						log((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i])) - cumsum_coal_rate[i] - 
						log(age_end - age_begin);
					assert(!std::isnan(num_e));
					num[e] = logsumexp(num[e], num_e);

					double x_begin = t_begin/age_end, x_end = t_end/age_end;
					double term1   = (x_begin*(age_end-t_begin)/inv_coal_rate_e + 1.0 - 2.0*(x_begin + inv_coal_rate_e/age_end));
					double term2   = (-x_end*(age_end -t_end)/inv_coal_rate_e   - 1.0 + 2.0*(x_end   + inv_coal_rate_e/age_end));
					double tmp     = term1;
					tmp           += exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) * term2;
					assert(!std::isnan(tmp));

					if(tmp < 0.0){
						denom[e] = log(0.0);
					}else{
						tmp        = log(tmp);
						tmp       += log(age_end) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
						assert(!std::isnan(tmp));
						denom_e    = log(DAF[i]) + tmp - log(age_end-age_begin);
						assert(!std::isnan(denom_e));
						denom_e    = logminusexp(denom_e, log(t_begin) + num_e) + log(num_lins[i]);
						assert(!std::isnan(denom_e));
						denom[e]   = logsumexp(denom[e], logsumexp(denom_e, log(cumsum_bl) + num_e));
						assert(!std::isnan(denom[e]));
					}

				}

			}else{
				break;
			}

			cumsum_bl += (t_end - t_begin)*num_lins[i];
			i++;
			if(i == num_age_bins-1) break;
		}

		if(normalising_constant == 1.0){
			normalising_constant = num[e];
		}else{
			normalising_constant = logsumexp(normalising_constant, num[e]);
		}

		if(i == num_age_bins-1 || DAF[i] == 0) break;

	}

	//std::cerr << std::endl;
	//e = 1;
	//std::cerr << exp(num[e]) << " " << exp(normalising_constant) << " " << exp(num[e] - normalising_constant) << std::endl;
	//std::cerr << exp(denom[e]) << " " << exp(normalising_constant) << " " << exp(denom[e] - normalising_constant) << std::endl;

	e = 0;
	if(!std::isinf(normalising_constant)){

		double integ = 1.0;
		i = 0;
		for(; e < std::min(num_epochs-1, ep_index[i] + 1); e++){

			double factor = 0.0;
			if(i < num_age_bins){
				while(ep_index[i] == e){
					factor += (t_int[i+1]-t_int[i]) * num_lins[i];
					i++;
					if(i == num_age_bins) break;
				}
			}

			num[e]   -= normalising_constant;
			denom[e] -= normalising_constant;
			num[e] = exp(num[e]);
			if(integ > 0.0){
				integ  -= num[e];
			}else{
				integ   = 0.0;
			}
			denom[e]  = exp(denom[e]);
			denom[e] += factor*integ;
			if(denom[e] < 0.0) denom[e] = 0.0;
			assert(!std::isnan(num[e]));
			assert(!std::isnan(denom[e]));
			assert(num[e] >= 0.0);
			assert(denom[e] >= 0.0);

		}
		if(ep_index[i] == num_epochs - 1){
			e = num_epochs-1;
			num[e]   -= normalising_constant;
			denom[e] -= normalising_constant;
			num[e]    = exp(num[e]);
			denom[e]  = exp(denom[e]);
			if(denom[e] < 0.0) denom[e] = 0.0;
			assert(!std::isnan(num[e]));
			assert(!std::isnan(denom[e]));
			assert(num[e] >= 0.0);
			assert(denom[e] >= 0.0);
			e++;
		}

	}
	for(; e < num_epochs; e++){
		num[e]   = exp(num[e]);
		denom[e] = exp(denom[e]);
	}

	return normalising_constant;

}

double
coal_EM_tree::EM_notshared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom){

	//prespecify t_int as a time grid, and ep_index according to t_int
	//input num_lin and DAF according to t_int
	//calculate integrals and sum up
	std::fill(num.begin(), num.end(), log(0.0));
	std::fill(denom.begin(), denom.end(), log(0.0));

	assert(age_begin <= age_end);
	bool age_identical = false;
	if(age_begin == age_end) age_identical = true;

	assert(num_lins.size() == num_age_bins);
	assert(DAF.size() == num_age_bins);

	double t_begin, t_end, inv_coal_rate_e;
	int current_e = -1;
	double normalising_constant = 1.0;

	int i = 0, e = 0;
	for(e = 0; e < num_epochs; e++){

		double num_e, denom_e;
		double cumsum_bl = 0.0;

		while(ep_index[i] == e){

			inv_coal_rate_e = 1.0/(num_lins[i] * coal_rates[e]);
			t_begin    = t_int[i];
			t_end      = t_int[i+1];

			if(t_end <= age_begin && DAF[i] < 1.0){
				//constant regime
				num_e    = log(1.0-DAF[i]) + logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
				assert(!std::isnan(num_e));
				num[e]   = logsumexp(num[e], num_e);

				denom_e  = log(1.0-DAF[i]) + 
					log((t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i])) - cumsum_coal_rate[i];
				if(std::isnan(denom_e)){
					std::cerr << t_begin << " " << t_end << " " << age_begin << " " << age_end << " " << inv_coal_rate_e << " " << cumsum_coal_rate[i] << " " << cumsum_coal_rate[i+1] << std::endl;
				}
				assert(!std::isnan(denom_e));
				denom_e  = logminusexp(denom_e, log(t_begin) + num_e) + log(num_lins[i]);
				assert(!std::isnan(denom_e));
				denom[e] = logsumexp(denom[e], logsumexp(denom_e, log(cumsum_bl) + num_e));
				assert(!std::isnan(denom[e]));
			}else if(t_begin >= age_begin && t_end <= age_end && !age_identical){

				if(1){
					//linear regime
					num_e  = log(DAF[i]) +
						log( (t_begin - age_begin + inv_coal_rate_e) + (age_begin - t_end - inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i])) - cumsum_coal_rate[i] - 
						log(age_end - age_begin);
					assert(!std::isnan(num_e));

					double x_begin = t_begin/age_end, x_end = t_end/age_end, x_age_begin = age_begin/age_end;
					double term1   = (x_begin*(t_begin - age_begin)/inv_coal_rate_e + 2.0*(x_begin + inv_coal_rate_e/age_end) - x_age_begin);
					double term2   = (-x_end*(t_end - age_begin)/inv_coal_rate_e    - 2.0*(x_end   + inv_coal_rate_e/age_end) + x_age_begin);
					double tmp     = term1;
					tmp           += exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) * term2;
					assert(!std::isnan(tmp));
					if(tmp < 0.0){
						denom[e] = log(0.0);
					}else{
						tmp        = log(tmp);
						tmp       += log(age_end) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
						assert(!std::isnan(tmp));
						denom_e    = log(DAF[i]) + tmp - log(age_end - age_begin);
						assert(!std::isnan(denom_e));
					}
				}else{
					num_e = log(0.0);
					denom_e = log(0.0);
				}

				//need to add linear regime to this
				num_e    = logsumexp(num_e, log(1.0-DAF[i]) + logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]));
				assert(!std::isnan(num_e));
				num[e]   = logsumexp(num[e], num_e);

				denom_e  = logsumexp(denom_e, log(1.0-DAF[i]) + 
						log((t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i])) - cumsum_coal_rate[i]);
				assert(!std::isnan(denom_e));
				denom_e    = logminusexp(denom_e, log(t_begin) + num_e) + log(num_lins[i]);
				assert(!std::isnan(denom_e));
				denom[e]   = logsumexp(denom[e], logsumexp(denom_e, log(cumsum_bl) + num_e));
				assert(!std::isnan(denom[e]));

			}else if(DAF[i] == 0){
				//constant regime and can coalesce anywhere
				num_e    = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
				assert(!std::isnan(num_e));
				num[e]   = logsumexp(num[e], num_e);

				denom_e  = log((t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i])) - cumsum_coal_rate[i];
				if(!std::isnan(denom_e)){
					assert(!std::isnan(denom_e));
					denom_e  = logminusexp(denom_e, log(t_begin) + num_e) + log(num_lins[i]);
					assert(!std::isnan(denom_e));
					denom[e] = logsumexp(denom[e], logsumexp(denom_e, log(cumsum_bl) + num_e));
					assert(!std::isnan(denom[e]));
				}
			}

			cumsum_bl += (t_end - t_begin)*num_lins[i];
			i++;
			if(i == num_age_bins-1) break;
		}

		if(normalising_constant == 1.0){
			normalising_constant = num[e];
		}else{
			normalising_constant = logsumexp(normalising_constant, num[e]);
		}

		if(i == num_age_bins-1) break;

	}

	if(!std::isinf(normalising_constant)){
		double integ = 1.0;	
		i = 0;
		for(e = 0; e < num_epochs-1; e++){

			double factor = 0.0;
			if(i < num_age_bins){
				while(ep_index[i] == e){
					factor += (t_int[i+1]-t_int[i]) * num_lins[i];
					i++;
					if(i == num_age_bins) break;
				}
			}

			num[e]   -= normalising_constant;
			denom[e] -= normalising_constant;
			num[e] = exp(num[e]);
			if(integ > 0.0){
				integ -= num[e];
			}else{
				integ  = 0.0;
			}
			denom[e]  = exp(denom[e]);
			denom[e] += factor*integ;
			if(denom[e] < 0.0) denom[e] = 0.0;
			assert(!std::isnan(num[e]));
			assert(!std::isnan(denom[e]));
			assert(num[e] >= 0.0);
			assert(denom[e] >= 0.0);
		}
		e = num_epochs-1;
		num[e]   -= normalising_constant;
		denom[e] -= normalising_constant;
		num[e]    = exp(num[e]);
		denom[e]  = exp(denom[e]);
		if(denom[e] < 0.0) denom[e] = 0.0;
		assert(!std::isnan(num[e]));
		assert(!std::isnan(denom[e]));
		assert(num[e] >= 0.0);
		assert(denom[e] >= 0.0);
	}else{
		std::fill(num.begin(), num.end(), 0.0);
		std::fill(denom.begin(), denom.end(), 0.0);
	}

	return normalising_constant;

}


////////////////////////////////////

double
coal_EM_tree_fast::logsumexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return logb;
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	if(loga > logb){
		long double res = loga + log1p(exp(logb - loga));
		return(res);
	}else{
		long double res = logb + log1p(exp(loga - logb));
		return(res);
	}

}

double
coal_EM_tree_fast::logminusexp(double loga, double logb){

	if(std::isinf(loga) || std::isnan(loga)){
		if(std::isinf(logb) || std::isnan(logb)){
			return log(0.0);
		}else{
			return log(0.0); //assuming small value
		}
	}
	if(std::isinf(logb) || std::isnan(logb)){
		if(std::isinf(loga) || std::isnan(loga)){
			return log(0.0);
		}else{
			return loga;
		}
	}

	if(loga < logb){
		return(log(0));
	}
	return(loga + log1p(-exp(logb - loga)));
	//return(loga + log(-expm1(logb-loga)));
	//return(loga + log(1.0 - exp(logb - loga)));

}

void
coal_EM_tree_fast::UpdateEpochs(std::vector<double>& iepochs){

	epochs = iepochs;
	num_epochs = epochs.size();
	factor.resize(num_epochs);

	it_tint = t_int.begin();
	it_ep_index = ep_index.begin();
	int bin = 0;
	int e   = 0;
	*it_ep_index = 0;
	it_tint++;
	for(bin = 1; bin < num_age_bins; bin++){
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

void
coal_EM_tree_fast::UpdateCoal(std::vector<double>& coal){

	it_coal = coal.begin();
	it_sum_coal = sum_coal_rate_tmpl.begin();
	it_inv_coal = inv_coal_rate_tmpl.begin();
	it_ep_index = ep_index.begin();
	it_tint = t_int.begin();
	it_tint_next = std::next(it_tint, 1);

	int e = *it_ep_index;
	for(; it_inv_coal != inv_coal_rate_tmpl.end();){
		if(e != num_epochs-1){
			if(e < *it_ep_index){
				it_coal++;
				e = *it_ep_index;
			}
		}
		*it_inv_coal = 1.0/(*it_coal);
		it_inv_coal++;
		it_ep_index++;
	}

	it_inv_coal = inv_coal_rate_tmpl.begin();
	for(; it_sum_coal != std::prev(sum_coal_rate_tmpl.end(),1);){
		*it_sum_coal = (*it_tint_next - *it_tint)/(*it_inv_coal);
		it_sum_coal++;
		it_inv_coal++;
		it_tint++;
		it_tint_next++;
	}

	it_sum_coal = sum_coal_rate_tmpl.begin();
	it1_exp_sum_coal = exp_sum_coal_rate_tint.begin();
	long double exp_coal;
	for(; it_sum_coal != sum_coal_rate_tmpl.end(); it_sum_coal++){
		exp_coal = exp(-*it_sum_coal);
		it2_exp_sum_coal = (*it1_exp_sum_coal).begin();	
		it2_exp_sum_coal_prev = it2_exp_sum_coal;
		//int k = 1;
		//*it2_exp_sum_coal = exp(-*it_sum_coal);
		*it2_exp_sum_coal = exp_coal;
		it2_exp_sum_coal++;
		//k++;
		for(;it2_exp_sum_coal != (*it1_exp_sum_coal).end();){  
			//*it2_exp_sum_coal = exp(-*it_sum_coal*k);
			*it2_exp_sum_coal = *it2_exp_sum_coal_prev * exp_coal;
			it2_exp_sum_coal++;
			it2_exp_sum_coal_prev++;
			//k++;
		}
		it1_exp_sum_coal++;
	}

}

void
coal_EM_tree_fast::UpdateTree(std::vector<float>& num_lins){

	inv_coal_rate_tint[0] = inv_coal_rate_tmpl[0]/num_lins[0];
	sum_coal_rate_tint[0] = exp_sum_coal_rate_tint[0][num_lins[0]-1];
	cumsum_coal_rate[0]   = 1.0;
	//sum_coal_rate_tint[0] = sum_coal_rate_tmpl[0]*num_lins[0];
	//cumsum_coal_rate[0]   = 0.0;
	int e = 0;
	for(int i = 1; i < t_int.size(); i++){
		if(e < ep_index[i]){
			e = ep_index[i];
		}
		inv_coal_rate_tint[i]   = inv_coal_rate_tmpl[i]/num_lins[i];
		sum_coal_rate_tint[i]   = exp_sum_coal_rate_tint[i][num_lins[i]-1];
		cumsum_coal_rate[i]     = cumsum_coal_rate[i-1] * sum_coal_rate_tint[i-1];
		//std::cerr << sum_coal_rate_tint[i] << " " << cumsum_coal_rate[i] << std::endl;
		//sum_coal_rate_tint[i]   = sum_coal_rate_tmpl[i]*num_lins[i];
		//cumsum_coal_rate[i]     = cumsum_coal_rate[i-1] + sum_coal_rate_tint[i-1];
		//sum_coal_rate_tint[i-1] = exp(-sum_coal_rate_tint[i-1]);
	}
	//sum_coal_rate_tint[t_int.size()-1] = exp(-sum_coal_rate_tint[t_int.size()-1]);

	//for(int i = 0; i < t_int.size(); i++){
	//	cumsum_coal_rate[i] = exp(-cumsum_coal_rate[i]);
	//}

	it_sum_coal = sum_coal_rate_tint.begin();
	it_inv_coal = inv_coal_rate_tint.begin();
	it_cumsum_coal = cumsum_coal_rate.begin();
	it_ep_index = ep_index.begin();
	it_tint = t_int.begin();
	it_tint_next = std::next(it_tint, 1);
	it_num_lins = num_lins.begin();
	it_f = f.begin();
	it_tf = tf.begin();
	it_tf_prec = tf_prec.begin();

	double cumsum_bl = 0.0;
	e = 0;
	for(; it_ep_index != ep_index.end();){
		if(e < *it_ep_index){
			e = *it_ep_index;
			cumsum_bl = 0.0;
		}

		*it_f        = (1.0 - *it_sum_coal) * (*it_cumsum_coal);
		*it_tf       = ((*it_tint + *it_inv_coal) - (*it_tint_next + *it_inv_coal) * (*it_sum_coal)) * (*it_cumsum_coal);
		assert(!std::isnan(*it_tf));
		*it_tf_prec  = *it_tf - (*it_tint)*(*it_f);
		*it_tf_prec *= (*it_num_lins);
		if(*it_tf_prec < 0){
			*it_tf_prec = 0.0;
		}
		*it_tf_prec += cumsum_bl * (*it_f);
		assert(!std::isnan(*it_tf_prec));
		cumsum_bl   += (*it_tint_next - *it_tint)*(*it_num_lins);

		it_sum_coal++;
		it_inv_coal++;
		it_cumsum_coal++;
		it_ep_index++;
		it_tint++;
		it_tint_next++;
		it_f++;
		it_tf++;
		it_tf_prec++;
		it_num_lins++;
	}

	it_ep_index = ep_index.begin();
	it_num_lins = num_lins.begin();
	it_tint = t_int.begin();
	it_tint_next = std::next(it_tint, 1);
	for(int e = 0; e < num_epochs; e++){

		factor[e] = 0.0;
		while(*it_ep_index == e){
			factor[e] += ((*it_tint_next) - (*it_tint)) * (*it_num_lins);

			it_ep_index++;
			it_num_lins++;
			it_tint++;
			it_tint_next++;
			if(it_tint_next == t_int.end()) break;
		}

	}

}

double 
coal_EM_tree_fast::Logl_shared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF){

	assert(age_begin <= age_end);
	bool age_identical = false;
	if(age_begin == age_end) age_identical = true;

	double t_begin, t_end, inv_coal_rate_e, num_e;
	double logl = 0.0;

	for(int i = 0; i < t_int.size();){

		t_begin         = t_int[i];
		t_end           = t_int[i+1];

		if(t_end <= age_begin){
			//constant regime
			logl   += DAF[i] * f[i];
		}else if(t_begin >= age_begin && t_end <= age_end && !age_identical){

			inv_coal_rate_e = inv_coal_rate_tint[i];
			//linear regime
			num_e  = DAF[i] * ((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * sum_coal_rate_tint[i]) * cumsum_coal_rate[i];
			num_e /= (age_end - age_begin);
			logl  += num_e;	

		}else{
			break;
		}

		i++;
		if(i == num_age_bins-1 || DAF[i] == 0) break;
	}

	if(logl == 0.0){
		return(0.0);
	}else{
		return(log(logl));
	}

}

double 
coal_EM_tree_fast::Logl_notshared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF){

	assert(age_begin <= age_end);
	bool age_identical = false;
	if(age_begin == age_end) age_identical = true;

	double t_begin, t_end, inv_coal_rate_e, num_e, denom_e;
	double logl = 0.0;

	for(int i = 0; i < t_int.size();){

		t_begin         = t_int[i];
		t_end           = t_int[i+1];

		if(t_end <= age_begin && DAF[i] < 1.0){
			//constant regime
			logl   += (1.0 - DAF[i]) * f[i];
		}else if(t_begin >= age_begin && t_end <= age_end && !age_identical){

			inv_coal_rate_e = inv_coal_rate_tint[i];
			//linear regime
			num_e  = DAF[i] * ((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * sum_coal_rate_tint[i]) * cumsum_coal_rate[i];
			num_e /= (age_end - age_begin);
			num_e  = std::max(0.0, f[i] - num_e);
			assert(!std::isinf(num_e));
			logl += num_e;

		}else if(DAF[i] == 0.0){
			logl   += f[i];
		}

		i++;
		if(i == num_age_bins-1) break;

	}

	if(logl == 0.0){
    return(0.0);
	}else{
	  return(log(logl));
	}
}

double
coal_EM_tree_fast::EM_shared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom){

	//prespecify t_int as a time grid, and ep_index according to t_int
	//input num_lin and DAF according to t_int
	//calculate integrals and sum up
	std::fill(num.begin(), num.end(), 0.0);
	std::fill(denom.begin(), denom.end(), 0.0);

	assert(age_begin <= age_end);
	bool age_identical = false;
	if(age_begin == age_end) age_identical = true;

	double t_begin, t_end, inv_coal_rate_e, num_e, denom_e;
	double normalising_constant = 0.0;

	int i = 0, e = 0;
	for(e = 0; e < num_epochs; e++){

		double cumsum_bl = 0.0;

		while(ep_index[i] == e){

			t_begin         = t_int[i];
			t_end           = t_int[i+1];

			if(t_end <= age_begin){
				//constant regime
				num[e]   += DAF[i] * f[i];
				denom[e] += DAF[i] * tf_prec[i];
				//if(e == 0) std::cerr << log(DAF[i] * tf_prec[i]) << " " << DAF[i] << " " << tf_prec[i] << "; ";
				//if(e == 1) std::cerr << log(DAF[i] * f[i]) << " ";
			}else if(t_begin >= age_begin && t_end <= age_end && !age_identical){

				inv_coal_rate_e = inv_coal_rate_tint[i];
				if(1){

					//linear regime
					num_e  = DAF[i] * ((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * sum_coal_rate_tint[i]) * cumsum_coal_rate[i];
					num_e /= (age_end - age_begin);
					assert(!std::isnan(num_e));
					num[e] += num_e;

					double x_begin = t_begin/age_end, x_end = t_end/age_end;
					double term1   = (x_begin*(age_end-t_begin)/inv_coal_rate_e + 1.0 - 2.0*(x_begin + inv_coal_rate_e/age_end));
					double term2   = (-x_end*(age_end -t_end)/inv_coal_rate_e   - 1.0 + 2.0*(x_end   + inv_coal_rate_e/age_end));
					double tmp     = term1;
					tmp           += sum_coal_rate_tint[i] * term2;
					//assert(!std::isnan(tmp));

					if(tmp < 0.0){
						denom[e] = 0.0;
					}else{
						tmp       *= age_end * inv_coal_rate_e * cumsum_coal_rate[i];
						//assert(!std::isnan(tmp));
						denom_e    = DAF[i]*tmp/(age_end-age_begin);

						//assert(!std::isnan(denom_e));
						denom_e   -= t_begin*num_e;
						denom_e   *= num_lins[i];
						//assert(!std::isnan(denom_e));
						denom_e   += cumsum_bl*num_e;
						denom[e]  += denom_e;
						assert(!std::isnan(denom_e));
					}

				}

			}else{
				break;
			}

			cumsum_bl += (t_end - t_begin)*num_lins[i];
			i++;
			if(i == num_age_bins-1) break;
		}

		normalising_constant += num[e];

		if(i == num_age_bins-1 || DAF[i] == 0) break;

	}

	//std::cerr << std::endl;
	//e = 1;
	//std::cerr << num[e] << " " << normalising_constant << " " << num[e]/normalising_constant << std::endl;
	//std::cerr << denom[e] << " " << normalising_constant << " " << denom[e]/normalising_constant << std::endl;

	if(normalising_constant > 0.0){

		double integ = 1.0;
		e = 0;
		for(; e < std::min(num_epochs-1, ep_index[i] + 1); e++){

			num[e]   /= normalising_constant;
			denom[e] /= normalising_constant;
			if(integ > 0.0){
				integ  -= num[e];
			}else{
				integ   = 0.0;
			}
			denom[e] += factor[e]*integ;
			if(denom[e] < 0.0) denom[e] = 0.0;
			if(!(num[e] >= 0.0)){

				for(int e = 0; e < num_epochs; e++){
          std::cerr << num[e] << " ";
				}
				std::cerr << std::endl;
				for(int e = 0; e < num_epochs; e++){
					std::cerr << denom[e] << " ";
				}
				std::cerr << std::endl;
				for(int e = 0; e < num_epochs; e++){
					std::cerr << inv_coal_rate_tmpl[e] << " ";
				}
				std::cerr << std::endl;

			}
			assert(!std::isnan(num[e]));
			assert(!std::isnan(denom[e]));
			assert(num[e] >= 0.0);
			assert(denom[e] >= 0.0);

		}
		if(ep_index[i] == num_epochs - 1){
			e = num_epochs-1;

			num[e]   /= normalising_constant;
			denom[e] /= normalising_constant;
			if(denom[e] < 0.0) denom[e] = 0.0;
			assert(!std::isnan(num[e]));
			assert(!std::isnan(denom[e]));
			assert(num[e] >= 0.0);
			assert(denom[e] >= 0.0);
			e++;
		}

	}

	return log(normalising_constant);

}

double
coal_EM_tree_fast::EM_notshared(double age_begin, double age_end, std::vector<float>& num_lins, std::vector<float>& DAF, std::vector<double>& num, std::vector<double>& denom){

	//prespecify t_int as a time grid, and ep_index according to t_int
	//input num_lin and DAF according to t_int
	//calculate integrals and sum up
	std::fill(num.begin(), num.end(), 0.0);
	std::fill(denom.begin(), denom.end(), 0.0);

	assert(age_begin <= age_end);
	bool age_identical = false;
	if(age_begin == age_end) age_identical = true;

	double t_begin, t_end, inv_coal_rate_e, num_e, denom_e;
	double normalising_constant = 0.0;

	int i = 0, e = 0;
	for(e = 0; e < num_epochs; e++){

		double cumsum_bl = 0.0;

		while(ep_index[i] == e){

			t_begin         = t_int[i];
			t_end           = t_int[i+1];

			if(t_end <= age_begin && DAF[i] < 1.0){
				//constant regime
				num[e]   += (1.0 - DAF[i]) * f[i];
				denom[e] += (1.0 - DAF[i]) * tf_prec[i];
			}else if(t_begin >= age_begin && t_end <= age_end && !age_identical){

				inv_coal_rate_e = inv_coal_rate_tint[i];
				if(1){

					//linear regime
					num_e  = DAF[i] * ((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * sum_coal_rate_tint[i]) * cumsum_coal_rate[i];
					num_e /= (age_end - age_begin);
					num_e  = std::max(0.0, f[i] - num_e);
					assert(!std::isinf(num_e));
					num[e] += num_e;

					double x_begin = t_begin/age_end, x_end = t_end/age_end;
					double term1   = (x_begin*(age_end-t_begin)/inv_coal_rate_e + 1.0 - 2.0*(x_begin + inv_coal_rate_e/age_end));
					double term2   = (-x_end*(age_end -t_end)/inv_coal_rate_e   - 1.0 + 2.0*(x_end   + inv_coal_rate_e/age_end));
					double tmp     = term1;
					tmp           += sum_coal_rate_tint[i] * term2;
					//assert(!std::isinf(tmp));

					if(tmp < 0.0){
						denom[e] = tf[i];
					}else{
						tmp       *= age_end * inv_coal_rate_e * cumsum_coal_rate[i];
						//assert(!std::isinf(tmp));
						denom_e    = DAF[i]*tmp/(age_end-age_begin);
						denom_e    = tf[i] - denom_e;
					}

					//assert(!std::isinf(denom_e));
					denom_e   -= t_begin*num_e;
					if(denom_e > 0){
						denom_e   *= num_lins[i];
						//assert(!std::isinf(denom_e));
						denom_e   += cumsum_bl*num_e;
						denom[e]  += denom_e;
					}
					assert(!std::isnan(denom_e));

				}

			}else if(DAF[i] == 0.0){
				num[e]   += f[i];
				denom[e] += tf_prec[i];
			}

			cumsum_bl += (t_end - t_begin)*num_lins[i];
			i++;
			if(i == num_age_bins-1) break;
		}

		normalising_constant += num[e];
		if(i == num_age_bins-1) break;

	}

	if(normalising_constant > 0.0){

		double integ = 1.0;	
		i = 0;
		for(e = 0; e < num_epochs-1; e++){

			num[e]   /= normalising_constant;
			denom[e] /= normalising_constant;
			if(integ > 0.0){
				integ -= num[e];
			}else{
				integ  = 0.0;
			}
			denom[e] += factor[e]*integ;
			if(denom[e] < 0.0) denom[e] = 0.0;
			assert(!std::isnan(num[e]));
			assert(!std::isnan(denom[e]));
			assert(num[e] >= 0.0);
			assert(denom[e] >= 0.0);
		}
		e = num_epochs-1;

		num[e]   /= normalising_constant;
		denom[e] /= normalising_constant;
		if(denom[e] < 0.0) denom[e] = 0.0;
		assert(!std::isnan(num[e]));
		assert(!std::isnan(denom[e]));
		assert(num[e] >= 0.0);
		assert(denom[e] >= 0.0);

	}

	return log(normalising_constant);

}


////////////////////////////////////

double
coal_EM_simplified::logsumexp(double iloga, double ilogb){

	long double loga = iloga;
	long double logb = ilogb;

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
		long double res = loga + log1p(exp(logb - loga));
		return(res);
	}else{
		long double res = logb + log1p(exp(loga - logb));
		return(res);
	}

}

double
coal_EM_simplified::logminusexp(double loga, double logb){

	//long double loga = iloga;
	//long double logb = ilogb;

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

	if(loga < logb){
		return(log(0));
	}
	return(loga + log1p(-exp(logb - loga)));
	//return(loga + log(-expm1(logb-loga)));
	//return(loga + log(1.0 - exp(logb - loga)));

}

void
coal_EM_simplified::get_tint(double age, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin){

	int e = 0;
	int i = 0;
	for(e = 0; e < num_epochs; e++){
		if(age < epochs[e] && i_begin == -1){
			t_int[i] = age;
			ep_index[i] = e-1;
			i_begin = i;
			i++;
		}
		t_int[i] = epochs[e];
		ep_index[i] = e;
		i++;
	}
	if(i_begin == -1){
		t_int[i] = age;
		ep_index[i] = num_epochs-1;
		i_begin = i;
		i++;
	}

}

void 
coal_EM_simplified::get_ABC(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B){

	std::vector<double> cumsum_coal_rate(t_int.size(), 0.0);
	for(int i = 1; i < t_int.size(); i++){
		cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * (t_int[i] - t_int[i-1]);
	}

	double t_begin, t_end, coal_rate_e;
	int i = 0;
	for(i = 0; i < t_int.size()-1; i++){

		t_begin     = t_int[i];
		t_end       = t_int[i+1];
		coal_rate_e = coal_rates[ep_index[i]];

		assert(t_end >= t_begin);

		if(coal_rate_e > 0 && t_end != 0 && t_end - t_begin > 0){
			A[i] = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);

			if(coal_rate_e == 0){ 
				B[i] = log(t_end - t_begin) - cumsum_coal_rate[i];
				if(t_begin > 0){
					B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], B[i]);
				}
				B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
			}else{
				if(t_begin == 0){
					if(1){
						assert(cumsum_coal_rate[i] == 0);
						B[i] = A[i] - log(coal_rate_e);
						B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
					}else{
						B[i] = log(exp(A[i] + cumsum_coal_rate[i+1])/coal_rate_e - t_end) - cumsum_coal_rate[i+1];
					}
				}else{
					if(1){
						if(log(t_begin)-cumsum_coal_rate[i] > log(t_end)-cumsum_coal_rate[i+1]){
							B[i] = logminusexp(log(t_begin)-cumsum_coal_rate[i], log(t_end)-cumsum_coal_rate[i+1]);
							B[i] = logsumexp(B[i], A[i] - log(coal_rate_e));
						}else{
							B[i] = logminusexp(log(t_end)-cumsum_coal_rate[i+1], log(t_begin)-cumsum_coal_rate[i]);
							B[i] = logminusexp(A[i] - log(coal_rate_e), B[i]);
						}
					}else{
						B[i] = log(t_begin * exp(-cumsum_coal_rate[i] + cumsum_coal_rate[i+1]) - t_end + exp(A[i] + cumsum_coal_rate[i+1])/coal_rate_e) - cumsum_coal_rate[i+1];
						//B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
						//B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
					}
				}
			}

			//to deal with numerical inaccuracies:
			if(t_begin > 0){
				if(B[i] > log(t_end)+A[i] || B[i] < log(t_begin) + A[i]){
					B[i] = log((t_begin + t_end)/2.0) + A[i];
				}
			}else if(B[i] > log(t_end)+A[i]){
				if(B[i] > log(t_end)+A[i]){
					B[i] = log(t_end/2) + A[i];
				}
			}

		}else{

			A[i] = log(0.0); 
			B[i] = log(0.0);

		}

		assert(!std::isnan(A[i]));
		assert(!std::isnan(B[i]));

	}
	i = t_int.size()-1;
	t_begin = t_int[i];
	coal_rate_e = coal_rates[ep_index[i]];
	if(coal_rate_e > 0){
		A[i] = -cumsum_coal_rate[i];
		//B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
		B[i] = log(t_begin * exp(-cumsum_coal_rate[i]) + exp(A[i])/coal_rate_e);
	}else{
		A[i] = log(0.0); 
		B[i] = log(0.0);
	}

	assert(!std::isnan(A[i]));
	assert(!std::isnan(B[i]));

}

void 
coal_EM_simplified::get_ABC_lazy(std::vector<double>& t_int, int i_begin, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B){

	std::vector<long double> cumsum_coal_rate(t_int.size(), 0.0);
	for(int i = 1; i < t_int.size(); i++){
		cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * (t_int[i] - t_int[i-1]);
	}

	double t_begin, t_end, coal_rate_e;
	int i = 0;
	int current_e = -1;
	for(i = 0; i < t_int.size()-1; i++){

		int e = ep_index[i];

		//current_e is epoch of i-1
		if(i != i_begin && i != i_begin - 1){

			A[i] = A_ep[e];
			B[i] = B_ep[e];

		}else{

			//std::cerr << e << " " << epochs[e] << std::endl;

			t_begin     = t_int[i];
			t_end       = t_int[i+1];
			coal_rate_e = coal_rates[ep_index[i]];

			assert(t_end >= t_begin);
			if(coal_rate_e > 0 && (t_end != 0 && t_end - t_begin > 0)){
				A[i] = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);

				if(coal_rate_e == 0){
					assert(1);
					B[i] = log(t_end - t_begin) - cumsum_coal_rate[i];
					if(t_begin > 0){
						B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], B[i]);
					}
					B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
				}else{
					if(t_begin == 0){
						if(1){
							B[i] = A[i] - log(coal_rate_e);
							B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
						}else{
							//B[i] = log(exp(A[i])/coal_rate_e - t_end * exp(-cumsum_coal_rate[i+1]));
							B[i] = log(exp(A[i]+cumsum_coal_rate[i+1])/coal_rate_e - t_end) - cumsum_coal_rate[i+1];
						}
					}else{
						//B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
						//B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
						if(1){
							if(log(t_begin)-cumsum_coal_rate[i] > log(t_end)-cumsum_coal_rate[i+1]){
								B[i] = logminusexp(log(t_begin)-cumsum_coal_rate[i], log(t_end)-cumsum_coal_rate[i+1]);
								B[i] = logsumexp(B[i], A[i] - log(coal_rate_e));
							}else{
								B[i] = logminusexp(log(t_end)-cumsum_coal_rate[i+1], log(t_begin)-cumsum_coal_rate[i]);
								B[i] = logminusexp(A[i] - log(coal_rate_e), B[i]);
							}
						}else{
							//B[i] = log(t_begin * exp(-cumsum_coal_rate[i]) - t_end * exp(-cumsum_coal_rate[i+1]) + exp(A[i])/coal_rate_e);
							B[i] = log(t_begin * exp(-cumsum_coal_rate[i]+cumsum_coal_rate[i+1]) - t_end + exp(A[i]+cumsum_coal_rate[i+1])/coal_rate_e) - cumsum_coal_rate[i+1];

							//B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
							//B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
							//std::cerr << B[i] << " " << log(t_begin)+A[i] << " " << log(t_end)+A[i] << std::endl; 
						}
					}

				}

				//to deal with numerical inaccuracies:
				if(t_begin > 0){
					if(B[i] > log(t_end)+A[i] || B[i] < log(t_begin)+A[i]){
						B[i] = log((t_begin + t_end)/2.0)+A[i];
					}
				}else if(B[i] > log(t_end)+A[i]){
					if(B[i] > log(t_end)+A[i]){
						B[i] = log(t_end/2.0)+A[i];
					}
				}

			}else{

				A[i] = log(0.0); 
				B[i] = log(0.0);

			}

		}

		if(std::isnan(B[i])){
			std::cerr << A[i] << " " << B[i] << " " << std::endl;
		}
		assert(!std::isnan(A[i]));
		assert(!std::isnan(B[i]));
		current_e = ep_index[i];

	}

	i = t_int.size()-1;
	if(current_e == num_epochs-1){

		t_begin = t_int[i];
		coal_rate_e = coal_rates[ep_index[i]];
		if(coal_rate_e > 0){
			A[i] = -cumsum_coal_rate[i];
			//B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
			B[i] = log(t_begin * exp(-cumsum_coal_rate[i]) + exp(A[i])/coal_rate_e);
		}else{
			A[i] = log(0.0); 
			B[i] = log(0.0);
		}
	}else{
		A[i] = A_ep[num_epochs-1];
		B[i] = B_ep[num_epochs-1];
	}

	assert(!std::isnan(A[i]));
	assert(!std::isnan(B[i]));

}

double
coal_EM_simplified::EM_shared(double age, std::vector<double>& num, std::vector<double>& denom){

	//calculate A, B, C, then f, and tf
	int num_epochs = epochs.size();

	//make vector of time intervals
	std::vector<double> t_int(num_epochs + 1, 0.0);
	std::vector<int> ep_index(num_epochs + 1, 0);
	int i_begin = -1; 
	get_tint(age, t_int, ep_index, i_begin);
	assert(i_begin > 0);

	//calculate A, B, C
	std::vector<double> A(num_epochs+1, 0.0), B(num_epochs+1, 0.0), f(num_epochs, 0.0), tf(num_epochs, 0.0);
	//get_ABC(t_int, ep_index, A, B);
	get_ABC_lazy(t_int, i_begin, ep_index, A, B);

	//now given A, B, C, I can calculate num and denom for each epoch
	int e = 0;
	int i = 0;
	double normalising_constant = 1.0;
	for(; e < ep_index[i_begin]+1; e++){

		assert(ep_index[i] == e);
		if(i < i_begin){

			f[e]  = A[i];
			if(epochs[e] == 0){
				tf[e] = B[i];
			}else{
				tf[e] = logminusexp(B[i], log(epochs[e]) + A[i]);
			}

			i++;

		}

		if(0){

			assert(f[e] <= 0.0);
			if(normalising_constant == 1.0){
				normalising_constant = f[e];
				assert(normalising_constant <= 0.0);
			}else{
				double tmp = normalising_constant;
				normalising_constant = logsumexp(normalising_constant, f[e]);
				if(!(normalising_constant <= 0.0)) std::cerr << "const: " << normalising_constant << " " << tmp << " " << f[e] << std::endl;
				assert(normalising_constant <= 0.0);
			}
		}

		if(i == i_begin) break;

	}

	normalising_constant = 0.0;
	for(e = 0; e < ep_index[i_begin]+1; e++){
		normalising_constant += exp(f[e]);
	}
	normalising_constant = log(normalising_constant);

	//normalise density
	for(e = 0; e < ep_index[i_begin]+1; e++){
		f[e]  -= normalising_constant;
		tf[e] -= normalising_constant;
		assert(!std::isnan(f[e]));
		assert(!std::isnan(tf[e]));
		f[e]   = exp(f[e]);
		tf[e]  = exp(tf[e]);
		assert(f[e] <= 1.0);
	}

	//calculate num and denom
	double integ = 1.0;
	for(e = 0; e < std::min(num_epochs - 1, ep_index[i_begin] + 1); e++){
		if(e > ep_index[i_begin]){
			assert(f[e] == 0.0);
			assert(tf[e] == 0.0);
		}
		num[e]  += f[e];
		integ   -= f[e];
		if(integ < 0 || e == ep_index[i_begin]) integ = 0.0;
		assert(integ >= 0.0);
		denom[e] += tf[e] + (epochs[e+1] - epochs[e]) * integ;

		assert(denom[e] <= (epochs[e+1] - epochs[e]));
	}
	if(ep_index[i_begin] == num_epochs-1){
		e = num_epochs-1;
		num[e] += f[e];
		denom[e] += tf[e];
	}
	//e = num_epochs - 1;
	//num[e] += f[e];
	//denom[e] += tf[e];

	return normalising_constant;

}

double
coal_EM_simplified::EM_notshared(double age, std::vector<double>& num, std::vector<double>& denom){

	//calculate A, B, C, then f, and tf
	int num_epochs = epochs.size();

	//make vector of time intervals
	std::vector<double> t_int(num_epochs + 1, 0.0);
	std::vector<int> ep_index(num_epochs + 1, 0);
	int i_begin = -1; 
	get_tint(age, t_int, ep_index, i_begin);
	assert(i_begin > 0);

	//calculate A, B, C
	std::vector<double> A(num_epochs+1, 0.0), B(num_epochs+1, 0.0), f(num_epochs, 0.0), tf(num_epochs, 0.0);
	//get_ABC(t_int, ep_index, A, B);
	get_ABC_lazy(t_int, i_begin, ep_index, A, B);

	//now given A, B, I can calculate num and denom for each epoch
	int e = ep_index[i_begin];
	int i = i_begin;
	double normalising_constant = 1.0;
	for(; e < num_epochs; e++){

		assert(ep_index[i] == e);

		f[e] = A[i];
		if(epochs[e] == 0){
			tf[e] = B[i];
		}else{
			tf[e] = logminusexp(B[i], log(epochs[e]) + A[i]);
		}
		i++;

		if(normalising_constant == 1.0){
			normalising_constant = f[e];
		}else{
			assert(normalising_constant <= 0.0);
			normalising_constant = logsumexp(normalising_constant, f[e]);
		}

		assert(!std::isnan(f[e]));

	}

	assert(!std::isnan(normalising_constant));

	//normalise density
	double sum = 0.0;
	for(e = ep_index[i_begin]; e < num_epochs; e++){
		f[e]  -= normalising_constant;
		tf[e] -= normalising_constant;
		f[e]   = exp(f[e]);
		tf[e]  = exp(tf[e]);
		sum += f[e];
		if(!(f[e] <= 1)) std::cerr << f[e] << " " << normalising_constant << std::endl;
		assert(f[e] <= 1.0);
	}
	//assert(std::fabs(sum - 1.0) < 1e-5);

	//calculate num and denom
	double integ = 1.0;
	for(e = 0; e < ep_index[i_begin]; e++){
		denom[e] += (epochs[e+1] - epochs[e]);
	}
	for(e = ep_index[i_begin]; e < num_epochs-1; e++){
		num[e]  += f[e];
		integ   -= f[e];
		if(integ < 0) integ = 0.0;
		assert(integ >= 0.0);
		denom[e] += tf[e] + (epochs[e+1] - epochs[e]) * integ;
		assert(denom[e] <= (epochs[e+1] - epochs[e]));
	}
	e = num_epochs - 1;
	num[e] += f[e];
	denom[e] += tf[e];

	return normalising_constant;

}

double
coal_EM_simplified::EM_shared_exact(double age, std::vector<double>& num, std::vector<double>& denom, double coal_rate){

	double lambda = coal_rate;	
	double norm     = logminusexp(0.0, -lambda*age);
	for(int e = 0; e < num_epochs-1; e++){
		double ep1 = epochs[e], ep2 = epochs[e+1];

		double A = std::min(ep1, age), B = std::min(age, ep2);
		double th_num   = logminusexp(-lambda*A, -lambda*B);


		double th_denom = logsumexp(log(A)-lambda*A, th_num - log(lambda));
		th_denom = logminusexp(th_denom, log(B) - lambda*B);
		th_denom = logminusexp(th_denom, log(ep1) + th_num);


		//double th_denom = log(A*exp(-lambda*A) + exp(th_num)/lambda - B*exp(-lambda*B) - ep1*exp(th_num));

		double rest = logminusexp(-lambda*B, -lambda*age);
		th_denom = logsumexp(th_denom, log(ep2 - ep1) + rest);

		th_num   -= norm;
		th_denom -= norm;

		th_num = exp(th_num);
		th_denom = exp(th_denom);
		num[e] = th_num;
		denom[e] = th_denom;
	}

	return(norm);

}

double
coal_EM_simplified::EM_notshared_exact(double age, std::vector<double>& num, std::vector<double>& denom, double coal_rate){

	double lambda = coal_rate;	
	double norm     = -lambda*age;
	for(int e = 0; e < num_epochs-1; e++){

		double ep1 = epochs[e], ep2 = epochs[e+1];
		double A = std::max(ep1, age), B = std::max(age, ep2);
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

		num[e] = th_num;
		denom[e] = th_denom;

	}
	int e = num_epochs-1;
	double ep1 = epochs[e];
	double A = std::max(ep1, age);
	double th_num = -lambda*A;
	double th_denom = log(A*exp(-lambda*A) + exp(th_num)/lambda - ep1*exp(th_num));
	th_num   -= norm;
	th_denom -= norm;
	num[e] = th_num;
	denom[e] = th_denom;

	return norm;

}




