#include "aDNA_EM.hpp"

double
logsumexp(double loga, double logb){

  if(std::isinf(loga) || std::isnan(loga)){
    if(std::isinf(logb) || std::isnan(logb)){
      return 0.0;
    }else{
      return logb;
    }
  }
  if(std::isinf(logb) || std::isnan(logb)){
    if(std::isinf(loga) || std::isnan(loga)){
      return 0.0;
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
logsumexp(double loga, double logb, double logc){

  double res = 0;
  if(loga > logb){
    if(logc > loga){
      //logc > loga > logb
      res = logc + log(1.0 + exp(loga - logc));
      res = res  + log(1.0 + exp(res - logb));
    }else if(logc > logb){
      //loga > logc > logb
      res = loga + log(1.0 + exp(logc - loga));
      res = res  + log(1.0 + exp(logb - res));
    }else{
      //loga > logb > logc
      res = loga + log(1.0 + exp(logb - loga));
      res = res  + log(1.0 + exp(res - logc));
    }
  }else{
    //logb > loga
    if(logc > logb){
      //logc > logb > loga
      res = logc + log(1.0 + exp(logb - logc));
      res = res  + log(1.0 + exp(res - loga));
    }else if(logc > loga){
      //logb > logc > loga
      res = logb + log(1.0 + exp(logc - logb));
      res = res  + log(1.0 + exp(loga - res));
    }else{
      //logb > loga > logc
      res = logb + log(1.0 + exp(loga - logb));
      res = res  + log(1.0 + exp(res - logc));
    }

  }

  return(res);

}

double
logminusexp(double loga, double logb){

  if(std::isinf(loga) || std::isnan(loga)){
    if(std::isinf(logb) || std::isnan(logb)){
      return 0.0;
    }else{
      return 0.0; //assuming small value
    }
  }
  if(std::isinf(logb) || std::isnan(logb)){
    if(std::isinf(loga) || std::isnan(loga)){
      return 0.0;
    }else{
      return loga;
    }
  }

  //assert(loga > logb);
  if(loga <= logb){
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
    cumsum_coal_rate[i] = cumsum_coal_rate[i-1] + coal_rates[ep_index[i]] * (t_int[i] - t_int[i-1]);
  }

  double t_begin, t_end, coal_rate_e;
  int i = 0;
  for(i = 0; i < num_epochs+2; i++){

    t_begin     = t_int[i];
    t_end       = t_int[i+1];
    coal_rate_e = coal_rates[ep_index[i]];

    //std::cerr << t_begin << " " << t_end << std::endl;
    assert(t_end >= t_begin);

    if(t_end != 0 && t_end - t_begin > 0){
      A[i] = -cumsum_coal_rate[i] + log(1-exp(cumsum_coal_rate[i]-cumsum_coal_rate[i+1]));

      if(t_begin == 0){
        B[i] = A[i] - log(coal_rate_e);
      }else{
        B[i] = logsumexp(log(t_begin) - cumsum_coal_rate[i], A[i] - log(coal_rate_e));
      }
      B[i] = logminusexp(B[i], log(t_end) - cumsum_coal_rate[i+1]);

      if(t_begin == 0){
        C[i] = log(2.0) - log(coal_rate_e) + B[i];
      }else{
        C[i] = logsumexp(2.0 * log(t_begin) - cumsum_coal_rate[i], log(2.0) - log(coal_rate_e) + B[i]);   
      }
      C[i] = logminusexp(C[i], 2.0 * log(t_end) - cumsum_coal_rate[i+1]);

    }else{

      A[i] = 0; //log0
      B[i] = 0; //log0
      C[i] = 0; //log0

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

  if(std::isinf(A[0]))  std::cerr << A[0] << std::endl;

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
        if(age_begin == 0){
          f[e] = logsumexp(f[e], B[i] - log(age_end));
        }else{
          f[e] = logsumexp(f[e], logminusexp(B[i], log(age_begin)+A[i]) - log(age_end - age_begin));
        }
        double tmp = 0;
        if(epochs[e] == 0){
          if(age_begin == 0){
            tmp = C[i];          
          }else{
            tmp = logminusexp(C[i], log(age_begin) + B[i]);
          }
        }else{
          assert(age_begin > 0);
          tmp   = logsumexp(C[i], log(age_begin) + log(epochs[e]) + A[i]);
          tmp   = logminusexp(tmp, log(epochs[e]+age_begin) + B[i]);
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
      //in linearly increasing regime
      if(age_begin == 0){
        f[e]  = B[i] - log(age_end);
        tf[e] = C[i] - log(age_end);
      }else{
        f[e]  = logminusexp(B[i], log(age_begin)+A[i]) - log(age_end - age_begin);
        if(epochs[e] == 0){
          tf[e] = logminusexp(C[i], log(age_begin) + B[i]) - log(age_end - age_begin);
        }else{
          tf[e] = logminusexp(logsumexp(C[i], log(age_begin) + log(epochs[e]) + A[i]), log(epochs[e]+age_begin) + B[i]) - log(age_end - age_begin);
        }
      }
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
    f[e]   = exp(f[e]);
    tf[e]  = exp(tf[e]);
  }

  //calculate num and denom
  std::fill(num.begin(), num.end(), 0.0);
  std::fill(denom.begin(), denom.end(), 0.0);
  for(e = 0; e < ep_index[i_tmrca]+1; e++){
    num[e] = f[e];
    denom[e] = tf[e];
    for(int e2 = e+1; e2 < num_epochs; e2++){
      denom[e] += (epochs[e+1] - epochs[e]) * f[e2];
    }
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
      //in P(D=1|t,coal) linear regime

      f[e]  = logminusexp(log(age_end) + A[i], B[i]) - log(age_end - age_begin);
      //std::cerr << i << " " << t_int[i] << " " << t_int[i+1] << " " << epochs[e] << " " << A[i] << " " << log(age_end) + A[i] << " " << B[i] << std::endl;
      if(epochs[e] == 0){
        tf[e] = logminusexp( log(age_end) + B[i], C[i] ) - log(age_end - age_begin);
      }else{
        //(age_end-t)*(t - epochs[e]) = -epochs[e]*age_end + (epochs[e] + age_end) t - t^2
        //std::cerr << log(epochs[e] + age_end) + B[i] << " " << log(epochs[e]) + log(age_end) + A[i] << std::endl;
        //std::cerr << logminusexp(log(epochs[e] + age_end) + B[i], log(epochs[e]) + log(age_end) + A[i]) << " " << C[i] << std::endl;
        //std::cerr << (logminusexp(log(epochs[e] + age_end) + B[i], log(epochs[e]) + log(age_end) + A[i]) - C[i]) << std::endl;
        tf[e] = logminusexp(logminusexp(log(epochs[e] + age_end) + B[i], log(epochs[e]) + log(age_end) + A[i]), C[i]) - log(age_end - age_begin);
      }
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
  }
  
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
  for(e = 0; e < ep_index[i_tmrca]+1; e++){
    num[e] = f[e];
    denom[e] = tf[e];
    for(int e2 = e+1; e2 < num_epochs; e2++){
      denom[e] += (epochs[e+1] - epochs[e]) * f[e2];
    }
  }

  return 0;

}



