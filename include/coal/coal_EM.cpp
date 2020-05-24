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
      if(coal_rates[e] > 0){
        num_e   = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
        denom_e = log( (t_begin + inv_coal_rate_e)/inv_coal_rate_e - (t_end + inv_coal_rate_e)/inv_coal_rate_e*exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) ) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
      }else{
        num_e = log_0;
        denom_e = log_0;
      }
      i++;
    }

    if(e >= ep_index[i_begin] && !times_identical){
      t_begin = t_int[i];
      t_end   = t_int[i+1];

      if(coal_rates[e] > 0){
        num[e] = log((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i])) - cumsum_coal_rate[i] - log(age_end - age_begin);

        double x_begin = t_begin/age_end, x_end = t_end/age_end;
        double term1 = (x_begin*(age_end-t_begin)/inv_coal_rate_e + 1.0 - 2.0*(x_begin + inv_coal_rate_e/age_end));
        double term2 = (-x_end*(age_end -t_end)/inv_coal_rate_e   - 1.0 + 2.0*(x_end   + inv_coal_rate_e/age_end));
        double tmp = term1;
        tmp       += exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) * term2;

        if(tmp < 0.0){
          denom[e] = log_0;
        }else{
          tmp        = log(tmp);
          tmp       += log(age_end) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
          denom[e]   = tmp - log(age_end-age_begin);
        }
      }else{
        num[e] = log_0;
        denom[e] = log_0;
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

  if(!std::isinf(normalising_constant) && !std::isnan(normalising_constant)){
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
  }else{
    normalising_constant = 0.0;
    std::fill(num.begin(), num.end(), 0.0);
    std::fill(denom.begin(), denom.end(), 0.0);
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
      if(coal_rate_e > 0){
        num[e]   = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
        denom[e] = log( (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) ) - cumsum_coal_rate[i];
        normalising_constant = num[e];
      }else{
        num[e] = log_0;
        denom[e] = log_0;
        normalising_constant = log_0;
      }
      i++;
      e++;
      for(; e < num_epochs; e++){
        num[e]   = A_ep[e];
        denom[e] = B_ep[e];
        normalising_constant = logsumexp(normalising_constant, num[e]);
      }
    }else{
      t_begin   = t_int[i];
      assert(coal_rate_e > 0);
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

      if(coal_rate_e > 0.0){
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
      }else{
        num[e]   = log(0.0);
        denom[e] = log(0.0);
      }

      i++;
      if(ep_index[i_end] == e){

        assert(i == i_end);
        if(e != num_epochs-1){
          t_begin  = t_int[i];
          t_end    = t_int[i+1];
          if(coal_rate_e > 0){
            num[e]   = logsumexp(num[e], logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]));
            denom[e] = logsumexp(denom[e], log( (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i]) ) - cumsum_coal_rate[i]);
          }else{
            num[e] = log_0;
            denom[e] = log_0;
          }
          i++;
        }else{
          t_begin   = t_int[i];
          if(coal_rate_e > 0){
            num[e]    = logsumexp(num[e], -cumsum_coal_rate[i]);
            denom[e]  = logsumexp(denom[e], log(t_begin + inv_coal_rate_e)-cumsum_coal_rate[i]);
          }else{
            num[e] = log_0;
            denom[e] = log_0;
          }
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

  if(!std::isinf(normalising_constant) && !std::isnan(normalising_constant)){
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
  }else{
    normalising_constant = 0.0;
    std::fill(num.begin(), num.end(), 0.0);
    std::fill(denom.begin(), denom.end(), 0.0);
  }
  return normalising_constant;

}

double
coal_EM::EM_shared_sampled(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom){

  std::uniform_real_distribution<double> dist_unif(0,1);   

  int num_epochs = epochs.size();
  std::fill(num.begin(), num.end(), log(0.0));
  std::fill(denom.begin(), denom.end(), log(0.0));
  double true_age_begin = age_begin, true_age_end = age_end;
  bool times_identical = true;
  bool age_identical = true;
  double nconst = log(0.0);

  int max_iter = 1;
  for(int iter = 0; iter < max_iter; iter++){

    age_begin = dist_unif(rng) * (true_age_end - true_age_begin) + true_age_begin;
    age_end   = age_begin;
    std::vector<double> num_s(num.size(), log(0.0)), denom_s(denom.size(), log(0.0));

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
        num_s[e]   = A_ep[i];
        denom_s[e] = B_ep[i];
        i++;

      }else if(e == ep_index[i_begin]){
        assert(i+1 == i_begin);
        t_begin  = t_int[i];
        t_end    = t_int[i+1];
        num_e   = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
        denom_e = log( (t_begin + inv_coal_rate_e)/inv_coal_rate_e - (t_end + inv_coal_rate_e)/inv_coal_rate_e*exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) ) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
        i++;
      }

      if(e >= ep_index[i_begin] && !times_identical){
        t_begin = t_int[i];
        t_end   = t_int[i+1];

        num_s[e] = log((age_end - t_begin - inv_coal_rate_e) + (t_end - age_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i])) - cumsum_coal_rate[i] - log(age_end - age_begin);

        double x_begin = t_begin/age_end, x_end = t_end/age_end;
        double term1 = (x_begin*(age_end-t_begin)/inv_coal_rate_e + 1.0 - 2.0*(x_begin + inv_coal_rate_e/age_end));
        double term2 = (-x_end*(age_end -t_end)/inv_coal_rate_e   - 1.0 + 2.0*(x_end   + inv_coal_rate_e/age_end));
        double tmp = term1;
        tmp       += exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) * term2;

        if(tmp < 0.0){
          denom_s[e] = log(0.0);
        }else{
          tmp        = log(tmp);
          tmp       += log(age_end) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
          denom_s[e]   = tmp - log(age_end-age_begin);
        }

        i++;
        while(ep_index[i] == e){
          i++;
          if(i == ep_index.size()) break;
        }
      }

      if(e == ep_index[i_begin]){
        if(times_identical){
          num_s[e]   = num_e;
          denom_s[e] = denom_e;
        }else{
          num_s[e]   = logsumexp(num_s[e], num_e);
          denom_s[e] = logsumexp(denom_s[e], denom_e);
        }
      }

      if(normalising_constant == 1.0){
        normalising_constant = num_s[e];
      }else{
        normalising_constant = logsumexp(normalising_constant, num_s[e]);
      }

      if(i == ep_index.size()) break;
    }

    for(int e = 0; e < num_epochs; e++){
      num[e]   = logsumexp(num[e], num_s[e]);
      denom[e] = logsumexp(denom[e], denom_s[e]);
    }
    nconst = logsumexp(nconst, normalising_constant);

  }

  double integ = 1.0;
  int e = 0;
  for(e = 0; e < num_epochs-1; e++){
    num[e]   -= nconst;
    denom[e] -= nconst;
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
  num[e]   -= nconst;
  denom[e] -= nconst;
  num[e]    = exp(num[e]);
  denom[e]  = exp(denom[e]);
  denom[e] -= epochs[e]*num[e];
  if(denom[e] < 0.0) denom[e] = 0.0;

  nconst -= log(max_iter);
  return nconst;

}

double
coal_EM::EM_notshared_sampled(double age_begin, double age_end, std::vector<double>& num, std::vector<double>& denom){

  std::uniform_real_distribution<double> dist_unif(0,1);   

  int num_epochs = epochs.size();
  std::fill(num.begin(), num.end(), log(0));
  std::fill(denom.begin(), denom.end(), log(0));
  double true_age_begin = age_begin, true_age_end = age_end;
  bool times_identical = true;
  bool age_identical = true;
  double nconst = log(0.0);

  int max_iter = 1;
  for(int iter = 0; iter < max_iter; iter++){

    age_begin = dist_unif(rng) * (true_age_end - true_age_begin) + true_age_begin;
    age_end   = age_begin;
    std::vector<double> num_s(num.size(), log(0.0)), denom_s(denom.size(), log(0.0));

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
        num_s[e]   = logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]);
        denom_s[e] = log( (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) ) - cumsum_coal_rate[i];
        normalising_constant = num_s[e];
        i++;
        e++;
        for(; e < num_epochs; e++){
          num_s[e]   = A_ep[e];
          denom_s[e] = B_ep[e];
          normalising_constant = logsumexp(normalising_constant, num_s[e]);
        }
      }else{
        t_begin   = t_int[i];
        num_s[e]    = -cumsum_coal_rate[i];
        denom_s[e]  = log(t_begin + inv_coal_rate_e)-cumsum_coal_rate[i];
        normalising_constant = num_s[e];
        i++;
      }

    }else{

      for(; e < ep_index[i_end]+1; e++){

        coal_rate_e = coal_rates[e];
        inv_coal_rate_e = 1.0/coal_rates[e];
        assert(ep_index[i] == e);
        t_begin = t_int[i];
        t_end   = t_int[i+1];

        num_s[e] = log( (t_begin - age_begin + inv_coal_rate_e) + (age_begin - t_end - inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i])) - cumsum_coal_rate[i] - log(age_end - age_begin);

        double x_begin = t_begin/age_end, x_end = t_end/age_end, x_age_begin = age_begin/age_end;

        double term1 = (x_begin*(t_begin - age_begin)/inv_coal_rate_e + 2.0*(x_begin + inv_coal_rate_e/age_end) - x_age_begin);
        double term2 = (-x_end*(t_end - age_begin)/inv_coal_rate_e    - 2.0*(x_end   + inv_coal_rate_e/age_end) + x_age_begin);

        double tmp = term1;
        tmp       += exp(-cumsum_coal_rate[i+1] + cumsum_coal_rate[i]) * term2;

        if(tmp < 0.0){
          denom_s[e] = log(0.0);
        }else{
          tmp        = log(tmp);
          tmp       += log(age_end) + log(inv_coal_rate_e) - cumsum_coal_rate[i];
          denom_s[e]   = tmp - log(age_end - age_begin);
        }
        i++;
        if(ep_index[i_end] == e){

          assert(i == i_end);
          if(e != num_epochs-1){
            t_begin  = t_int[i];
            t_end    = t_int[i+1];
            num_s[e]   = logsumexp(num_s[e], logminusexp(-cumsum_coal_rate[i], -cumsum_coal_rate[i+1]));
            denom_s[e] = logsumexp(denom_s[e], log( (t_begin + inv_coal_rate_e) - (t_end + inv_coal_rate_e) * exp(-cumsum_coal_rate[i+1]+cumsum_coal_rate[i]) ) - cumsum_coal_rate[i]);
            i++;
          }else{
            t_begin   = t_int[i];
            num_s[e]    = logsumexp(num_s[e], -cumsum_coal_rate[i]);
            denom_s[e]  = logsumexp(denom_s[e], log(t_begin + inv_coal_rate_e)-cumsum_coal_rate[i]);
            i++;
          }

        }

        if(normalising_constant == 1.0){
          normalising_constant = num_s[e];
        }else{
          normalising_constant = logsumexp(normalising_constant, num_s[e]);
        }
      }

      for(; e < num_epochs; e++){
        num_s[e] = A_ep[e];
        denom_s[e] = B_ep[e];
        normalising_constant = logsumexp(normalising_constant, num_s[e]);
      }

    }

    for(int e = 0; e < num_epochs; e++){
      num[e]   = logsumexp(num[e], num_s[e]);
      denom[e] = logsumexp(denom[e], denom_s[e]);
    }
    nconst = logsumexp(nconst, normalising_constant);

  }

  double integ = 1.0;
  int e = 0;
  for(e = 0; e < num_epochs-1; e++){
    num[e]   -= nconst;
    denom[e] -= nconst;
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
  num[e]   -= nconst;
  denom[e] -= nconst;
  num[e]    = exp(num[e]);
  denom[e]  = exp(denom[e]);
  denom[e] -= epochs[e]*num[e];
  if(denom[e] < 0.0) denom[e] = 0.0;

  nconst -= log(max_iter);
  return nconst;

}


