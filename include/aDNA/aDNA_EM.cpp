#include "aDNA_EM.hpp"

double
aDNA_EM::logsumexp(double loga, double logb){

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
aDNA_EM::logminusexp(double loga, double logb){

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

void
aDNA_EM::get_tint(float age_begin, float age_end, float tmrca, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin, int& i_end, int& i_tmrca){

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
    if(tmrca     < epochs[e] && i_tmrca == -1){
      t_int[i] = tmrca;
      ep_index[i] = e-1;
      i_tmrca = i;
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
  if(i_tmrca == -1){
    t_int[i] = tmrca;
    ep_index[i] = num_epochs-1;
    i_tmrca = i;
    i++;
  }

}

void 
aDNA_EM::get_ABC(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B, std::vector<double>& C){

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
          B[i] = A[i] - log(coal_rate_e);
        }else{
          B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
        }
        B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
      }

			//to deal with numerical inaccuracies:
			if(t_begin > 0){
				if(B[i] > log(t_end)+A[i] || B[i] < log(t_begin) + A[i]){
					//B[i] = 0.5 * (log(t_end) + A[i] + log(t_begin) + A[i]);
          //std::cerr << log(t_end) - log(t_begin) << " " << B[i] - log(t_end) - A[i] << " " << log(t_begin) + A[i] - B[i] << std::endl;
          //std::cerr << "correcting1" << std::endl;
          B[i] = log((t_begin + t_end)/2.0) + A[i];
				}
			}else if(B[i] > log(t_end)+A[i]){
				if(B[i] > log(t_end)+A[i]){
          //std::cerr << "correcting2" << std::endl;
					B[i] = log(t_end) + A[i];
				}
			}

      if(coal_rate_e == 0){
        C[i] = log(t_end*t_end-t_begin*t_begin) - cumsum_coal_rate[i];
        if(t_begin > 0){
          C[i] = logsumexp(2.0 * log(t_begin) - cumsum_coal_rate[i], C[i]);
        }
        C[i] = logminusexp(C[i], 2.0 * log(t_end) - cumsum_coal_rate[i+1]);
      }else{
        if(t_begin == 0){
          C[i] = log(2.0) - log(coal_rate_e) + B[i];
        }else{
          C[i] = logsumexp(2.0 * log(t_begin) - cumsum_coal_rate[i], log(2.0) - log(coal_rate_e) + B[i]);   
        }
        C[i] = logminusexp(C[i], 2.0 * log(t_end) - cumsum_coal_rate[i+1]);
      }
			//to deal with numerical inaccuracies:
			if(t_begin > 0){
				if(C[i] > log(t_end)+B[i] || C[i] < log(t_begin)+B[i]){
					//C[i] = 0.5 * (log(t_end)+B[i]+log(t_begin)+B[i]);
          //std::cerr << "correcting3" << std::endl;
					C[i] = log((t_begin + t_end)/2.0)+B[i];
				}
			}else if(C[i] > log(t_end)+B[i]){
        //std::cerr << "correcting4" << std::endl;
				C[i] = log(t_end)+B[i];
			}
    }else{

      A[i] = log(0.0); 
      B[i] = log(0.0);
      C[i] = log(0.0);

    }

    assert(!std::isnan(A[i]));
    assert(!std::isnan(B[i]));
    assert(!std::isnan(C[i]));

  }
  i = t_int.size()-1;
  t_begin = t_int[i];
  coal_rate_e = coal_rates[ep_index[i]];
  if(coal_rate_e > 0){
    A[i] = -cumsum_coal_rate[i];
    B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
    C[i] = logsumexp(2.0*log(t_begin) - cumsum_coal_rate[i], log(2) - log(coal_rate_e) + B[i]);
  }else{
    A[i] = log(0.0); 
    B[i] = log(0.0);
    C[i] = log(0.0);
  }
  assert(!std::isnan(A[i]));
  assert(!std::isnan(B[i]));
  assert(!std::isnan(C[i]));

  if(std::isnan(A[0]))  std::cerr << A[0] << std::endl;

}

void 
aDNA_EM::get_ABC_lazy(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B, std::vector<double>& C){

  std::vector<double> cumsum_coal_rate(t_int.size(), 0.0);
  for(int i = 1; i < t_int.size(); i++){
    cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * (t_int[i] - t_int[i-1]);
  }

  double t_begin, t_end, coal_rate_e;
  int i = 0;
  int current_e = -1;
  for(i = 0; i < t_int.size()-1; i++){

		//std::cerr << "--- " << i << std::endl;

    int e = ep_index[i];
    if(e == ep_index[i+1]-1 && e == current_e+1){

      A[i] = A_ep[e];
      B[i] = B_ep[e];
      C[i] = C_ep[e];

    }else{


      t_begin     = t_int[i];
      t_end       = t_int[i+1];
      coal_rate_e = coal_rates[ep_index[i]];

      assert(t_end >= t_begin);
      if(coal_rate_e > 0 && (t_end != 0 && t_end - t_begin > 0)){
        A[i] = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);

        if(coal_rate_e == 0){ 
          B[i] = log(t_end - t_begin) - cumsum_coal_rate[i];
          if(t_begin > 0){
            B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], B[i]);
          }
          B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
        }else{
          if(t_begin == 0){
            B[i] = A[i] - log(coal_rate_e);
          }else{
            B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
          }
					B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]); 
        }

				//to deal with numerical inaccuracies:
				if(t_begin > 0){
					if(B[i] > log(t_end)+A[i] || B[i] < log(t_begin)+A[i]){
            //std::cerr << "correcting5" << std::endl;
            //std::cerr << log(t_begin)+A[i] << " " << B[i] << " " << log(t_end)+A[i] << std::endl;
						//B[i] = 0.5 * (log(t_end)+A[i]+log(t_begin)+A[i]);
            B[i] = log((t_begin + t_end)/2.0)+A[i];
					}
				}else if(B[i] > log(t_end)+A[i]){
					if(B[i] > log(t_end)+A[i]){
            //std::cerr << "correcting6" << std::endl;
            //std::cerr << B[i] << " " << log(t_end)+A[i] << std::endl;
						B[i] = log(t_end/2.0)+A[i];
					}
				}

        if(coal_rate_e == 0){
          C[i] = log(t_end*t_end-t_begin*t_begin) - cumsum_coal_rate[i];
          if(t_begin > 0){
            C[i] = logsumexp(2.0 * log(t_begin) - cumsum_coal_rate[i], C[i]);
          }
          C[i] = logminusexp(C[i], 2.0 * log(t_end) - cumsum_coal_rate[i+1]);
        }else{
          if(t_begin == 0){
            C[i] = logminusexp(log(2.0) - log(coal_rate_e) + B[i], 2.0 * log(t_end) - cumsum_coal_rate[i+1]);
          }else{
            if(0){
              C[i] = logsumexp(2.0 * log(t_begin/t_end) - cumsum_coal_rate[i], log(2.0) - log(coal_rate_e) + B[i] - 2.0*log(t_end));
              C[i] = logminusexp(C[i],-cumsum_coal_rate[i+1]);
              C[i] += 2.0*log(t_end);
            }else{
              C[i] = logsumexp(2.0 * log(t_begin) - cumsum_coal_rate[i], log(2.0) - log(coal_rate_e) + B[i]);
              C[i] = logminusexp(C[i], 2.0*log(t_end)- cumsum_coal_rate[i+1]);
            }
          }
				}
				//to deal with numerical inaccuracies:
				if(t_begin > 0){
					if(C[i] > log(t_end)+B[i] || C[i] < log(t_begin)+B[i]){
            //std::cerr << "correcting7" << std::endl;
            //std::cerr << log(t_begin)+B[i] << " " << C[i] << " " << log(t_end)+B[i] << std::endl;
						//C[i] = 0.5 * (log(t_end)+B[i] + log(t_begin)+B[i]);
            C[i] = log((t_begin + t_end)/2.0)+B[i];
					}
				}else if(C[i] > log(t_end)+B[i]){
          //std::cerr << "correcting8" << std::endl;
          C[i] = log(t_end/2.0) + B[i];
				}

				if(0){
				std::cerr << "check:" << std::endl;
			  std::cerr << log(t_begin) + A[i] << " " << B[i] << " " << log(t_end) + A[i] << std::endl;
				std::cerr << log(t_begin) + B[i] << " " << C[i] << " " << log(t_end) + B[i] << std::endl;
				std::cerr << 2.0*log(t_begin)+A[i] << " " << C[i] << " " << 2.0*log(t_end) + A[i] << std::endl;
				assert(B[i] <= log(t_end)   + A[i]);
				assert(B[i] >= log(t_begin) + A[i]);
				assert(C[i] <= log(t_end)   + B[i]);
				assert(C[i] >= log(t_begin) + B[i]);

				std::cerr << "------------" << std::endl; 
				}

      }else{

        A[i] = log(0.0); 
        B[i] = log(0.0);
        C[i] = log(0.0);

      }

    }

    assert(!std::isnan(A[i]));
    assert(!std::isnan(B[i]));
    assert(!std::isnan(C[i]));

    current_e = ep_index[i];

  }

  i = t_int.size()-1;
  if(current_e == num_epochs-1){

    t_begin = t_int[i];
    coal_rate_e = coal_rates[ep_index[i]];
    if(coal_rate_e > 0){
      A[i] = -cumsum_coal_rate[i];
      B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
      C[i] = logsumexp(2.0*log(t_begin) - cumsum_coal_rate[i], log(2) - log(coal_rate_e) + B[i]);
    }else{
      A[i] = log(0.0); 
      B[i] = log(0.0);
      C[i] = log(0.0);
    }
  }else{
    A[i] = A_ep[num_epochs-1];
    B[i] = B_ep[num_epochs-1];
    C[i] = C_ep[num_epochs-1];
  }

  assert(!std::isnan(A[i]));
  assert(!std::isnan(B[i]));
  assert(!std::isnan(C[i]));

  if(std::isnan(A[0]))  std::cerr << A[0] << std::endl;

}

double
aDNA_EM::EM_shared(float age_begin, float age_end, float tmrca, std::vector<float>& num, std::vector<float>& denom){

  //calculate A, B, C, then f, and tf
  std::fill(num.begin(), num.end(), 0);
  std::fill(denom.begin(), denom.end(), 0);

  int num_epochs = epochs.size();

  assert(age_begin < age_end);
  if(tmrca < age_end){
    tmrca = age_end;
  }
  assert(age_end <= tmrca);

  //make vector of time intervals
  std::vector<double> t_int(num_epochs + 3, 0.0);
  std::vector<int> ep_index(num_epochs+3, 0);
  int i_begin = -1, i_end = -1, i_tmrca = -1; 
  get_tint(age_begin, age_end, tmrca, t_int, ep_index, i_begin, i_end, i_tmrca);
  assert(ep_index[i_end] == ep_index[i_end - 1]);

  //calculate A, B, C
  std::vector<double> A(num_epochs+3, 0.0), B(num_epochs+3, 0.0), C(num_epochs+3, 0.0), f(num_epochs, 0.0), tf(num_epochs, 0.0);
  //get_ABC(t_int, ep_index, A, B, C);
  get_ABC_lazy(t_int, ep_index, A, B, C);

  //now given A, B, C, I can calculate num and denom for each epoch
  int e = 0;
  int i = 0;
  double normalising_constant = 0.0;
  for(; e < ep_index[i_end]+1; e++){

    assert(ep_index[i] == e);
    if(i < i_begin){

      //in P(D=1|t,coal) = 1 regime
      if(ep_index[i_begin] == e){
      
        //this will end in this epoch
        f[e] = A[i];
        if(epochs[e] == 0){
          tf[e] = B[i];
        }else{
          tf[e] = logminusexp(B[i], log(epochs[e]) + A[i]);
        }
        i++;
        //linear regime
        f[e] = logsumexp(f[e], logminusexp(log(age_end) + A[i], B[i]) - log(age_end - age_begin));

        double tmp = 0;
        tmp   = logminusexp(log(age_end + epochs[e]) + B[i], C[i]);
        if(epochs[e] > 0){
          tmp   = logminusexp(tmp, log(age_end) + log(epochs[e]) + A[i]);
        }
        tf[e] = logsumexp(tf[e], tmp - log(age_end - age_begin));
        i++;
        while(ep_index[i] == e) i++;

      }else{
     
        //this will not end in this epoch
        f[e]  = A[i];
        if(epochs[e] == 0){
          tf[e] = B[i];
        }else{
          tf[e] = logminusexp(B[i], log(epochs[e]) + A[i]);
        }
        i++;
     
      }

    }else if(i < i_end){
      
      //in linearly decreasing regime
      f[e]   = logminusexp(log(age_end) + A[i], B[i]) - log(age_end - age_begin);
      tf[e]  = logminusexp(log(age_end + epochs[e]) + B[i], C[i]);
      if(epochs[e] > 0){
        tf[e]  = logminusexp(tf[e], log(age_end) + log(epochs[e]) + A[i]);
      }
      tf[e] -= log(age_end - age_begin);
      i++;
      while(ep_index[i] == e) i++;

    }

    if(normalising_constant == 0.0){
      normalising_constant = f[e];
    }else{
      normalising_constant = logsumexp(normalising_constant, f[e]);
    }

  }

  double sum = 0.0, diff = 1.0;
  //normalise density
  for(e = 0; e < ep_index[i_end]+1; e++){
    f[e]  -= normalising_constant;
    tf[e] -= normalising_constant;
    assert(!std::isnan(f[e]));
    assert(!std::isnan(tf[e]));
    f[e]   = exp(f[e]);
    tf[e]  = exp(tf[e]);

    sum  += f[e];
    diff -= f[e];
    assert(f[e] <= 1.0);
  }

  //std::cerr << sum << " " << (sum <= 1.0) << " " << (sum == 1.0) << std::endl;
  //std::cerr << diff << " " << (diff >= 0.0) << " " << (diff == 0.0) << std::endl;

  //calculate num and denom
  std::fill(num.begin(), num.end(), 0.0);
  std::fill(denom.begin(), denom.end(), 0.0);
  double integ = 1.0;
  for(e = 0; e < num_epochs-1; e++){
		if(e > ep_index[i_end]){
      assert(f[e] == 0.0);
			assert(tf[e] == 0.0);
	  }
    num[e]   = f[e];
    integ   -= f[e];
    //if(integ < 0) std::cerr << e << " " << num_epochs-1 << " " << integ << " " << f[e] << " " << f[e+1] << std::endl;
    if(integ < 0) integ = 0.0;
    assert(integ >= 0.0);
    denom[e] = tf[e] + (epochs[e+1] - epochs[e]) * integ;
  }

  return normalising_constant;

}

double
aDNA_EM::EM_notshared(float age_begin, float age_end, float tmrca, std::vector<float>& num, std::vector<float>& denom){

  //calculate A, B, C, then f, and tf
  std::fill(num.begin(), num.end(), 0);
  std::fill(denom.begin(), denom.end(), 0);
  int num_epochs = epochs.size();

  assert(age_begin < age_end);
  if(tmrca < age_end) tmrca = age_end;
  assert(age_end <= tmrca);
  bool tmrca_eq_end = false;
  if(tmrca == age_end) tmrca_eq_end = true;

  //make vector of time intervals
  std::vector<double> t_int(num_epochs + 3, 0.0);
  std::vector<int> ep_index(num_epochs+3, 0);
  int i_begin = -1, i_end = -1, i_tmrca = -1; 
  get_tint(age_begin, age_end, tmrca, t_int, ep_index, i_begin, i_end, i_tmrca);

  //calculate A, B, C
  std::vector<double> A(num_epochs+3, 0.0), B(num_epochs+3, 0.0), C(num_epochs+3, 0.0), f(num_epochs, 0.0), tf(num_epochs, 0.0);
  //get_ABC(t_int, ep_index, A, B, C);
  get_ABC_lazy(t_int, ep_index, A, B, C);

  //now given A, B, C, I can calculate num and denom for each epoch
  int e = ep_index[i_begin];
  int i = i_begin;
  double normalising_constant = 0.0;
  for(; e < ep_index[i_tmrca]+1; e++){

    assert(ep_index[i] == e);

    if(i < i_end){
      //in P(D=0|t,coal) linear regime

      if(age_begin == 0){
        f[e]  = B[i] - log(age_end);
      }else{
        f[e]  = logminusexp(B[i], log(age_begin) + A[i]) - log(age_end - age_begin);
      }
      if(epochs[e] > 0 || age_begin > 0){
        if(epochs[e] > 0 && age_begin > 0){
          tf[e] = logsumexp(C[i], log(age_begin) + log(epochs[e]) + A[i]);
        }else{
          tf[e] = C[i];
        }
        tf[e] = logminusexp(tf[e], log(age_begin + epochs[e]) + B[i]);
      }else{
        tf[e] = C[i];
      }
      tf[e] -= log(age_end - age_begin);
      i++;

      if(ep_index[i_end] == e){
        //this will end in this epoch
        assert(i == i_end);
        if(!tmrca_eq_end){
          f[e] = logsumexp(f[e], A[i]);
          if(epochs[e] == 0){
            tf[e] = logsumexp(tf[e], B[i]);
          }else{
            tf[e] = logsumexp(tf[e], logminusexp(B[i], log(epochs[e]) + A[i]));
          }
        }
        i++;
      }

    }else if(i < i_tmrca){
      //in P(D=1|t,coal) = 1 regime

      f[e] = A[i];
      if(epochs[e] == 0){
        tf[e] = B[i];
      }else{
        tf[e] = logminusexp(B[i], log(epochs[e]) + A[i]);
      }
      i++;

    }

    if(normalising_constant == 0.0){
      normalising_constant = f[e];
    }else{
      normalising_constant = logsumexp(normalising_constant, f[e]);
    }

    if(e == ep_index[i_tmrca]) assert(i == i_tmrca);

  }

  float exp_t = 0.0;
  //normalise density
  for(e = ep_index[i_begin]; e < ep_index[i_tmrca]+1; e++){
    f[e]  -= normalising_constant;
    tf[e] -= normalising_constant;
    f[e]   = exp(f[e]);
    tf[e]  = exp(tf[e]);
		assert(f[e] <= 1.0);
  }
  
  //calculate num and denom
  std::fill(num.begin(), num.end(), 0.0);
  std::fill(denom.begin(), denom.end(), 0.0);
  double integ = 1.0;
	for(e = 0; e < num_epochs-1; e++){
		if(e > ep_index[i_tmrca] || e < ep_index[i_begin]){
			assert(f[e] == 0.0);
			assert(tf[e] == 0.0);
		}
    num[e]   = f[e];
    integ   -= f[e];
    if(integ < 0) integ = 0.0;
    assert(integ >= 0.0);
    denom[e] = tf[e] + (epochs[e+1] - epochs[e]) * integ;
  }

  return normalising_constant;

}

////////////////////////////////////

double
aDNA_EM_simplified::logsumexp(double iloga, double ilogb){

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
		//if(res > 0.0 && loga < 0.0 && logb < 0.0) return loga;
		return(res);
		//return(loga);
  }else{
		long double res = logb + log1p(exp(loga - logb));
		//if(res > 0.0 && loga < 0.0 && logb < 0.0) return logb;
		return(res);
		//return(logb);
  }

}

double
aDNA_EM_simplified::logminusexp(double loga, double logb){

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

  //assert(loga > logb);
  if(loga < logb){
    //std::cerr << loga << " " << logb << std::endl;
    return(log(0));
  }
  return(loga + log1p(-exp(logb - loga)));
  //return(loga + log(-expm1(logb-loga)));
  //return(loga + log(1.0 - exp(logb - loga)));

}

void
aDNA_EM_simplified::get_tint(double age, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin){

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
aDNA_EM_simplified::get_ABC(std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B){

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
					//std::cerr << "debug" << std::endl;
          B[i] = log((t_begin + t_end)/2.0) + A[i];
				}
			}else if(B[i] > log(t_end)+A[i]){
				if(B[i] > log(t_end)+A[i]){
					//std::cerr << "debug" << std::endl;
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
		B[i] = log(t_begin * exp(cumsum_coal_rate[i]) + exp(A[i])/coal_rate_e);
  }else{
    A[i] = log(0.0); 
    B[i] = log(0.0);
  }

	if(std::isnan(B[i])){
		std::cerr << A[i] << " " << B[i] << " " << std::endl;
	}
  assert(!std::isnan(A[i]));
  assert(!std::isnan(B[i]));

}

void 
aDNA_EM_simplified::get_ABC_lazy(std::vector<double>& t_int, int i_begin, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B){

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
    //if(e == ep_index[i+1]-1 && e == current_e+1){
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
						//std::cerr << B[i] << " " << log(t_end)+A[i] << " " << log(t_begin)+A[i] << std::endl;
						//std::cerr << "debug" << std::endl;
            B[i] = log((t_begin + t_end)/2.0)+A[i];
					}
				}else if(B[i] > log(t_end)+A[i]){
					if(B[i] > log(t_end)+A[i]){
						//std::cerr << "debug" << std::endl;
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
			B[i] = log(t_begin * exp(cumsum_coal_rate[i]) + exp(A[i])/coal_rate_e);
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

  if(std::isnan(A[0]))  std::cerr << A[0] << std::endl;

}

double
aDNA_EM_simplified::EM_shared(double age, std::vector<double>& num, std::vector<double>& denom){

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
				//tf[e] = logminusexp(B[i], log(epochs[e]) + A[i]);
				tf[e] = log(exp(B[i]-A[i]) - epochs[e]) + A[i];
			}

			i++;

    }

		if(0){
    
		assert(f[e] <= 0.0);
		if(normalising_constant == 1.0){
      normalising_constant = f[e];
			assert(normalising_constant <= 0.0);
    }else{
			//if(!(normalising_constant <= 0.0)) std::cerr << "const: " << normalising_constant << std::endl;
      //assert(normalising_constant <= 0.0);
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
  //std::fill(num.begin(), num.end(), 0.0);
  //std::fill(denom.begin(), denom.end(), 0.0);
  double integ = 1.0;
  //for(e = 0; e < num_epochs-1; e++){
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
  //e = num_epochs - 1;
  //num[e] += f[e];
  //denom[e] += tf[e];

  return normalising_constant;

}

double
aDNA_EM_simplified::EM_notshared(double age, std::vector<double>& num, std::vector<double>& denom){

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
      //tf[e] = logminusexp(B[i], log(epochs[e]) + A[i]);
			tf[e] = log(exp(B[i] - A[i]) - epochs[e]) + A[i];
    }
    i++;

    if(normalising_constant == 1.0){
      normalising_constant = f[e];
    }else{
      assert(normalising_constant <= 0.0);
			//std::cerr << normalising_constant << " " << f[e] << " " << age << " " << i_begin << " " << ep_index[i_begin] << " " << e << " " << num_epochs << " " << epochs[e] << " " << A_ep[e] << " " << A[i-1] << " " << A[i] << " " << A[i_begin] << std::endl;
			//std::cerr << coal_rates[e] << std::endl;
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
aDNA_EM_simplified::EM_shared_exact(double age, std::vector<double>& num, std::vector<double>& denom){

  double lambda = 1/40000.0;	
	double norm     = logminusexp(0.0, -lambda*age);
	for(int e = 0; e < num_epochs-1; e++){
		double ep1 = epochs[e], ep2 = epochs[e+1];

		double A = std::min(ep1, age), B = std::min(age, ep2);
		double th_num   = logminusexp(-lambda*A, -lambda*B);

		/*
		double th_denom = logsumexp(log(A)-lambda*A, th_num - log(lambda));
		th_denom = logminusexp(th_denom, log(B) - lambda*B);
		th_denom = logminusexp(th_denom, log(ep1) + th_num);
		*/

		double th_denom = log(A*exp(-lambda*A) + exp(th_num)/lambda - B*exp(-lambda*B) - ep1*exp(th_num));

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
aDNA_EM_simplified::EM_notshared_exact(double age, std::vector<double>& num, std::vector<double>& denom){

	double lambda = 1/40000.0;	
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




