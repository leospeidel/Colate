#include "catch.hpp"

#include "data.hpp"
#include "anc.hpp"
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


TEST_CASE("test EM expectation step"){

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

  double initial_coal_rate = 1.0/30000.0;
  Data data(1, 1);
  std::vector<std::vector<double>> coal_rates(data.N), coal_rates_num(data.N), coal_rates_denom(data.N);
  for(int i = 0; i < data.N; i++){
    coal_rates[i].resize(num_epochs);
    coal_rates_num[i].resize(num_epochs);
    coal_rates_denom[i].resize(num_epochs);
    std::fill(coal_rates[i].begin(), coal_rates[i].end(), initial_coal_rate);
  }

  aDNA_EM EM(epochs, coal_rates[0]);

  int i = 0;

  int bp_ref, bp_input;
  std::vector<char> sequence_ref(data.N), sequence_input(data.N);
  std::vector<float> num(num_epochs), denom(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1

  float age_begin = 8e4, age_end = 9e4, tmrca;
  //REQUIRE(EM.EM_shared(    age_begin, age_end, 10000000, num, denom) == 0);
  //REQUIRE(EM.EM_notshared( age_begin, age_end, 10000000, num, denom) == 0);

  age_begin = 1e6/28;
  age_end   = 5e7/28;
  tmrca     = 1e8/28;
  //std::cerr << age_begin << " " << age_end << std::endl;
  REQUIRE(EM.EM_shared(    age_begin, age_end, tmrca, num, denom) == 0);
  REQUIRE(!std::isnan(denom[0]));
  REQUIRE(EM.EM_notshared( age_begin, age_end, tmrca, num, denom) == 0);
  REQUIRE(!std::isnan(denom[0]));


  double x = 0, y = 1, lambda = 1/30000.0;
  int e = 20;
  x = 0; 
  y = 1;
  for(int i = 0; i < 6; i++){
  
    y *= 10;
    //y = std::max(10.0, epochs[i]-100);
    for(int e = num_epochs-2; e >= 0; e--){
      double ep1 = epochs[e], ep2 = epochs[e+1];
      double A = std::min(std::max(ep1, x), y), B = std::max(std::min(y, ep2), x), C = std::min(y, std::max(ep2,x));
      double norm     = logminusexp(-lambda*x, -lambda*y);
      double th_num   = logminusexp(-lambda*A, -lambda*B)-norm;
      
      double th_denom = logsumexp(log(A)-lambda*A,-log(lambda)-lambda*A);
      th_denom = logminusexp(th_denom, log(B)-lambda*B);
      th_denom = logminusexp(th_denom, -log(lambda)-lambda*B);
      th_denom = logminusexp(th_denom, log(ep1)+logminusexp(-lambda*A,-lambda*B));
      th_denom = logsumexp(th_denom, log(ep2-ep1)+logminusexp(-lambda*C,-lambda*y));
			th_denom -= norm;

      EM.EM_shared(    y, 1.0001*y, 1e9, num, denom);

			if(0){
      std::cerr << e << " " << x << " " << y << std::endl;
      std::cerr << A << " " << B << " " << C << " " << x << " " << y << std::endl;
			std::cerr << norm << std::endl;
      std::cerr << num[e] << " " << denom[e] << " " << denom[e]/num[e] << std::endl;
      std::cerr << exp(th_num) << " " << exp(th_denom) << " " << exp(th_denom-th_num) << std::endl;
			}

			th_num = exp(th_num);
			th_denom = exp(th_denom);

      if(!std::isnan(num[e]) && !std::isnan(th_num)){
        if(th_num != 0 && num[e] != 0){
          REQUIRE(std::fabs(num[e] - th_num)/th_num <= 0.01);
        }else{
          REQUIRE(std::fabs(num[e]) <= 0.01);
        }
      }
      if(!std::isnan(denom[e]) && !std::isnan(th_denom)){
        if(th_denom != 0 && denom[e] != 0){
          REQUIRE(std::fabs(denom[e] - th_denom)/th_denom <= 0.01);
        }else{
          REQUIRE(std::fabs(denom[e]) <= 0.01);
        }
      }

    }
  }

  y = 1e9;
  x = 1; 
  for(int i = 0; i < 6; i++){
  
    x *= 10;
    //x = std::max(10.0, epochs[i]-100);
    for(int e = num_epochs-2; e >= 0; e--){
      double ep1 = epochs[e], ep2 = epochs[e+1];
      double A = std::min(std::max(ep1, x), y), B = std::max(std::min(y, ep2), x), C = std::min(y, std::max(ep2,x));
      double norm     = logminusexp(-lambda*x, -lambda*y);
      double th_num   = logminusexp(-lambda*A, -lambda*B)-norm;
      
      //std::cerr << "v1: " << exp(logminusexp(logminusexp( -log(lambda), log(B)-lambda*B), -log(lambda)-lambda*B)-norm) << std::endl;
     
      double th_denom = logsumexp(log(A)-lambda*A,-log(lambda)-lambda*A);
      th_denom = logminusexp(th_denom, log(B)-lambda*B);
      th_denom = logminusexp(th_denom, -log(lambda)-lambda*B);
      th_denom = logminusexp(th_denom, log(ep1)+logminusexp(-lambda*A,-lambda*B));
      th_denom = logsumexp(th_denom, log(ep2-ep1)+logminusexp(-lambda*C,-lambda*y));
			th_denom -= norm;

      EM.EM_notshared(    x, 1.0001*x, 1e9, num, denom);

			if(0){
      std::cerr << e << " " << x << " " << y << std::endl;
      std::cerr << A << " " << B << " " << C << " " << x << " " << y << std::endl;
			std::cerr << "n " << norm << std::endl;
      std::cerr << num[e] << " " << denom[e] << " " << denom[e]/num[e] << std::endl;
      std::cerr << exp(th_num) << " " << exp(th_denom) << " " << exp(th_denom-th_num) << std::endl;
			}

			th_num = exp(th_num);
			th_denom = exp(th_denom);

      if(!std::isnan(num[e]) && !std::isnan(th_num)){
        if(th_num != 0 && num[e] != 0){
          REQUIRE(std::fabs(num[e] - th_num)/th_num <= 0.01);
        }else{
          REQUIRE(std::fabs(num[e]) <= 0.01);
        }
      }
      if(!std::isnan(denom[e]) && !std::isnan(th_denom)){
        if(th_denom != 0 && denom[e] != 0){
          REQUIRE(std::fabs(denom[e] - th_denom)/th_denom <= 0.01);
        }else{
          REQUIRE(std::fabs(denom[e]) <= 0.01);
        }
      }

    }
  }

  if(1){
  //need to work out some theoretical results here
  REQUIRE(EM.EM_shared(10, 100, 100, num, denom) == 0);
  REQUIRE(EM.EM_notshared(10, 100, 100, num, denom) == 0);

  REQUIRE(EM.EM_shared(10, 100, 10000, num, denom) == 0);
  REQUIRE(EM.EM_notshared(10, 100, 10000, num, denom) == 0);

  REQUIRE(EM.EM_shared(0, 100, 10000, num, denom) == 0);
  REQUIRE(EM.EM_notshared(0, 100, 10000, num, denom) == 0);

  REQUIRE(EM.EM_shared(1000, 5000, 10000, num, denom) == 0);
  REQUIRE(EM.EM_notshared(1000, 5000, 10000, num, denom) == 0);

  REQUIRE(EM.EM_shared(1000, 1001, 10000, num, denom) == 0);
  REQUIRE(EM.EM_notshared(1000, 1001, 10000, num, denom) == 0);

  REQUIRE(EM.EM_shared(epochs[3], 1001, 10000, num, denom) == 0);
  REQUIRE(EM.EM_notshared(epochs[3], 1001, 10000, num, denom) == 0);

  REQUIRE(EM.EM_shared(epochs[3], epochs[5], 10000, num, denom) == 0);
  REQUIRE(EM.EM_notshared(epochs[3], epochs[5], 10000, num, denom) == 0);

  aDNA_EM EM2(epochs, coal_rates[0]);

  for(int i = 0; i < data.N; i++){
    coal_rates[i][10] = 0;
  }
  REQUIRE(EM2.EM_shared(10, 100, 10000, num, denom) == 0);
  REQUIRE(EM2.EM_notshared(10, 100, 10000, num, denom) == 0);
  }

}

