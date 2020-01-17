#ifndef ADNA_EM_HPP
#define ADNA_EM_HPP

#include "gzstream.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"

#include <iostream>
#include <string>

int EM_shared(std::vector<double>& coal_rates, float age_begin, float age_end, float tmrca, std::vector<float>& epochs, std::vector<float>& num, std::vector<float>& denom);
int EM_notshared(std::vector<double>& coal_rates, float age_begin, float age_end, float tmrca, std::vector<float>& epochs, std::vector<float>& num, std::vector<float>& denom); 

#endif //ADNA_EM_HPP

