#include "aDNA_EM.hpp"

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

void
get_tint(std::vector<float>& epochs, int num_epochs, float age_begin, float age_end, float tmrca, std::vector<double>& t_int, std::vector<int>& ep_index, int& i_begin, int& i_end, int& i_tmrca){

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
    t_int[i] == age_begin;
    ep_index[i] = num_epochs-1;
    i_begin = i;
    i++;
  }
  if(i_end == -1){
    t_int[i] == age_end;
    ep_index[i] = num_epochs-1;
    i_end = i;
    i++;
  }
  if(i_tmrca == -1){
    t_int[i] == tmrca;
    ep_index[i] = num_epochs-1;
    i_tmrca = i;
    i++;
  }
  assert(i == num_epochs + 3);

}

void 
get_ABC(int num_epochs, std::vector<double>& coal_rates, std::vector<double>& t_int, std::vector<int>& ep_index, std::vector<double>& A, std::vector<double>& B, std::vector<double>& C){

  std::vector<double> cumsum_coal_rate(num_epochs+3, 0.0);
  for(int i = 1; i < num_epochs+3; i++){
    cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i-1]] * (t_int[i] - t_int[i-1]);
  }

  double t_begin, t_end, coal_rate_e;
  int i = 0;
  for(i = 0; i < num_epochs+2; i++){

    t_begin     = t_int[i];
    t_end       = t_int[i+1];
    coal_rate_e = coal_rates[ep_index[i]];


    //std::cerr << i << " " << t_begin << " " << t_end << " " << coal_rate_e << " " << cumsum_coal_rate[i] << " " << cumsum_coal_rate[i+1] << std::endl;
    assert(t_end >= t_begin);

    if(t_end != 0 && t_end - t_begin > 0){
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
    }else{

      A[i] = log(0.0); 
      B[i] = log(0.0);
      C[i] = log(0.0);

    }

    assert(!std::isnan(A[i]));
    assert(!std::isnan(B[i]));
    assert(!std::isnan(C[i]));

  }
  i = num_epochs+2;
  t_begin = t_int[i];
  coal_rate_e = coal_rates[ep_index[i]];
  A[i] = -cumsum_coal_rate[i];
  B[i] = log(t_begin) - cumsum_coal_rate[i] + A[i] - log(coal_rate_e);
  C[i] = 2.0*log(t_begin) - cumsum_coal_rate[i] + log(2) - log(coal_rate_e) + B[i];

  assert(!std::isnan(A[i]));
  assert(!std::isnan(B[i]));
  assert(!std::isnan(C[i]));

  if(std::isnan(A[0]))  std::cerr << A[0] << std::endl;

}

int
EM_shared(std::vector<double>& coal_rates, float age_begin, float age_end, float tmrca, std::vector<float>& epochs, std::vector<float>& num, std::vector<float>& denom){

  //calculate A, B, C, then f, and tf

  int num_epochs = epochs.size();

  assert(age_begin < age_end);
  if(tmrca < age_end) tmrca = age_end;
  assert(age_end <= tmrca);

  //make vector of time intervals
  std::vector<double> t_int(num_epochs + 3, 0.0);
  std::vector<int> ep_index(num_epochs+3, 0);
  int i_begin = -1, i_end = -1, i_tmrca = -1; 
  get_tint(epochs, num_epochs, age_begin, age_end, tmrca, t_int, ep_index, i_begin, i_end, i_tmrca);

  //calculate A, B, C
  std::vector<double> A(num_epochs+3, 0.0), B(num_epochs+3, 0.0), C(num_epochs+3, 0.0), f(num_epochs, 0.0), tf(num_epochs, 0.0);
  get_ABC(num_epochs, coal_rates, t_int, ep_index, A, B, C);

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

  //normalise density
  for(e = 0; e < ep_index[i_end]+1; e++){
    f[e]  -= normalising_constant;
    tf[e] -= normalising_constant;
    assert(!std::isnan(f[e]));
    assert(!std::isnan(tf[e]));
    f[e]   = exp(f[e]);
    tf[e]  = exp(tf[e]);
  }

  //calculate num and denom
  std::fill(num.begin(), num.end(), 0.0);
  std::fill(denom.begin(), denom.end(), 0.0);
  double integ = 1.0;
  for(e = 0; e < ep_index[i_tmrca]+1; e++){
    num[e]   = f[e];
    integ   -= f[e];
    if(integ < 0) integ = 0.0;
    denom[e] = tf[e] + (epochs[e+1] - epochs[e]) * integ;
    //std::cerr << epochs[e] << " " << num[e] << " " << denom[e] << std::endl;
  }
  //std::cerr << std::endl;
 
  return 0;

}

int
EM_notshared(std::vector<double>& coal_rates, float age_begin, float age_end, float tmrca, std::vector<float>& epochs, std::vector<float>& num, std::vector<float>& denom){

  //calculate A, B, C, then f, and tf

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
  get_tint(epochs, num_epochs, age_begin, age_end, tmrca, t_int, ep_index, i_begin, i_end, i_tmrca);

  //calculate A, B, C
  std::vector<double> A(num_epochs+3, 0.0), B(num_epochs+3, 0.0), C(num_epochs+3, 0.0), f(num_epochs, 0.0), tf(num_epochs, 0.0);
  get_ABC(num_epochs, coal_rates, t_int, ep_index, A, B, C);

  //now given A, B, C, I can calculate num and denom for each epoch
  int e = ep_index[i_begin];
  int i = i_begin;
  double normalising_constant = 0.0;
  for(; e < ep_index[i_tmrca]+1; e++){

    //std::cerr << e << " " << epochs[e] << " " << i << " " << i_end << std::endl;
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

  double tmp1 = 0, tmp2 = 0, tmp3 = 0;
  //normalise density
  for(e = ep_index[i_begin]; e < ep_index[i_tmrca]+1; e++){
    //std::cerr << epochs[e] << " " << f[e] << " " << tf[e] << std::endl;
    f[e]  -= normalising_constant;
    tf[e] -= normalising_constant;
    f[e]   = exp(f[e]);
    tf[e]  = exp(tf[e]);
    //std::cerr << epochs[e] << " " << f[e] << std::endl;
  }
  //std::cerr << std::endl;
  
  for(e = 0; e < num_epochs; e++){

    tmp1  += f[e];
    tmp2  += tf[e] + epochs[e]*f[e];

    tmp3 = 0;
    for(int e2 = e+1; e2 < num_epochs; e2++){
      tmp3 += f[e2];
    }

    //std::cerr << epochs[e] << "\t" << f[e] << "\t" << tmp3 << std::endl;// << " " << tf[e] + (epochs[e+1] - epochs[e])*tmp3 << std::endl;

  }

  //std::cerr << std::endl;
  //std::cerr << tmp1 << " " << tmp2 << std::endl;

  //calculate num and denom
  std::fill(num.begin(), num.end(), 0.0);
  std::fill(denom.begin(), denom.end(), 0.0);
  double integ = 1.0;
  for(e = 0; e < ep_index[i_tmrca]+1; e++){
    num[e]   = f[e];
    integ   -= f[e];
    if(integ < 0) integ = 0.0;
    denom[e] = tf[e] + (epochs[e+1] - epochs[e]) * integ;
  }

  return 0;

}



