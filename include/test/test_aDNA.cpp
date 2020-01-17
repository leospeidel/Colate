#include "catch.hpp"

#include "data.hpp"
#include "anc.hpp"
#include "aDNA_EM.hpp"

TEST_CASE("test EM expectation step"){

  //decide on epochs
  float years_per_gen = 28.0;
  int num_epochs = 30;
  num_epochs++;
  std::vector<float> epochs(num_epochs);
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

  int i = 0;

  int bp_ref, bp_input;
  std::vector<char> sequence_ref(data.N), sequence_input(data.N);
  std::vector<float> num(num_epochs), denom(num_epochs); //temporary variables storing numerator and demonmitor of MLE for SNP given D=0 or D=1

  //need to work out some theoretical results here
  REQUIRE(EM_shared(coal_rates[i], 10, 100, 100, epochs, num, denom) == 0);
  REQUIRE(EM_notshared(coal_rates[i], 10, 100, 100, epochs, num, denom) == 0);

  REQUIRE(EM_shared(coal_rates[i], 10, 100, 10000, epochs, num, denom) == 0);
  REQUIRE(EM_notshared(coal_rates[i], 10, 100, 10000, epochs, num, denom) == 0);

  REQUIRE(EM_shared(coal_rates[i], 0, 100, 10000, epochs, num, denom) == 0);
  REQUIRE(EM_notshared(coal_rates[i], 0, 100, 10000, epochs, num, denom) == 0);

  REQUIRE(EM_shared(coal_rates[i], 1000, 5000, 10000, epochs, num, denom) == 0);
  REQUIRE(EM_notshared(coal_rates[i], 1000, 5000, 10000, epochs, num, denom) == 0);

}

