#include <iostream>

#include "gzstream.hpp"
#include "data.hpp"
#include "sample.hpp"
#include "anc.hpp"
#include "mutations.hpp"
#include "cxxopts.hpp"

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

struct threepop_stats{

  std::string pop1, pop2, pop3;
  int a, b, c1, c2;
  double mean_ac1_bc2 = 0.0, sd_ac1_bc2 = 0.0;
  double mean_ac2_bc1 = 0.0, sd_ac2_bc1 = 0.0;
  double mean_ac1b_c2 = 0.0, sd_ac1b_c2 = 0.0;
  double mean_bc1a_c2 = 0.0, sd_bc1a_c2 = 0.0;
  double mean_ac2b_c1 = 0.0, sd_ac2b_c1 = 0.0;
  double mean_bc2a_c1 = 0.0, sd_bc2a_c1 = 0.0;

  double mean_ac1bc2 = 0.0, mean_ac2bc1 = 0.0;
  double sd_ac1bc2   = 0.0, sd_ac2bc1 = 0.0;

  double count_ac1_bc2 = 0.0;
  double count_ac2_bc1 = 0.0;
  double count_ac1b_c2 = 0.0;
  double count_bc1a_c2 = 0.0;
  double count_ac2b_c1 = 0.0;
  double count_bc2a_c1 = 0.0;

};

void      
get_C(Node& node, Sample& sample, int pop, std::vector<int>& C){

  if(node.child_left != NULL){
   
    std::vector<int> C_left, C_right; 
    get_C(*node.child_left, sample, pop, C_left);
    get_C(*node.child_right, sample, pop, C_right);

    C.resize(C_left.size() + C_right.size());
    std::vector<int>::iterator it_C_left = C_left.begin(), it_C_right = C_right.begin(), it_C = C.begin();

    for(; it_C != C.end();){
      
      if(it_C_left != C_left.end() && it_C_right != C_right.end()){
        if(*it_C_left < *it_C_right){
          *it_C = *it_C_left;
          it_C_left++;
        }else{
          *it_C = *it_C_right;
          it_C_right++;
        }
      }else if(it_C_left != C_left.end()){
        *it_C = *it_C_left;
        it_C_left++;
      }else if(it_C_right != C_right.end()){
        *it_C = *it_C_right;
        it_C_right++;
      }

      it_C++;
    }

  }else{
    if(sample.group_of_haplotype[node.label] == pop){
      C.push_back(node.label);
    }
  }


}

void
assign_ftree_3pop(int a, int b, std::vector<int>& C1, std::vector<int>& C2, std::vector<int>& C3, std::vector<threepop_stats>& ftree, std::vector<int>& convert_pop, int num_a, int num_b, int num_c){

  //C1 is clade of pop samples underneath a
  //C2 is clade of pop samples underneath b
  //C3 is clade of pop samples above (a,b)
  
  int entry_index = (convert_pop[a] * num_b + convert_pop[b]) * ((num_c*(num_c-1.0))/2.0);

  //iterate through all pairs of c's
  std::vector<int>::iterator it_C1 = C1.begin();
  std::vector<int>::iterator it_C2 = C2.begin();
  std::vector<int>::iterator it_C3 = C3.begin();

  for(it_C1 = C1.begin(); it_C1 != C1.end(); it_C1++){
    
    for(it_C2 = C2.begin(); it_C2 != C2.end(); it_C2++){
      //std::cerr << *it_C1 << " " << *it_C2 << std::endl;
      assert(*it_C1 != *it_C2);
      if(*it_C1 > *it_C2){
        ftree[entry_index + (convert_pop[*it_C1]*(convert_pop[*it_C1]-1.0))/2.0 + convert_pop[*it_C2]].count_ac1_bc2++;
      }else{
        ftree[entry_index + (convert_pop[*it_C2]*(convert_pop[*it_C2]-1.0))/2.0 + convert_pop[*it_C1]].count_ac2_bc1++;
      }
    }

    for(it_C3 = C3.begin(); it_C3 != C3.end(); it_C3++){     
      assert(*it_C1 != *it_C3);
      if(*it_C1 > *it_C3){
        ftree[entry_index + (convert_pop[*it_C1]*(convert_pop[*it_C1]-1.0))/2.0 + convert_pop[*it_C3]].count_ac1b_c2++;
        //ftree[(convert_pop[*it_C1]*(convert_pop[*it_C1]-1.0))/2.0 + convert_pop[*it_C3]].mean_ac2b_c1++;
      }else{
        ftree[entry_index + (convert_pop[*it_C3]*(convert_pop[*it_C3]-1.0))/2.0 + convert_pop[*it_C1]].count_ac2b_c1++;
        //ftree[(convert_pop[*it_C3]*(convert_pop[*it_C3]-1.0))/2.0 + convert_pop[*it_C1]].mean_ac2b_c1++;
      }
    }

  }

  for(it_C2 = C2.begin(); it_C2 != C2.end(); it_C2++){
    for(it_C3 = C3.begin(); it_C3 != C3.end(); it_C3++){
   
      assert(*it_C2 != *it_C3);
      if(*it_C2 > *it_C3){
        ftree[entry_index + (convert_pop[*it_C2]*(convert_pop[*it_C2]-1.0))/2.0 + convert_pop[*it_C3]].count_bc1a_c2++;
        //ftree[(convert_pop[*it_C2]*(convert_pop[*it_C2]-1.0))/2.0 + convert_pop[*it_C3]].mean_bc2a_c1++;
      }else{
        ftree[entry_index + (convert_pop[*it_C3]*(convert_pop[*it_C3]-1.0))/2.0 + convert_pop[*it_C2]].count_bc2a_c1++;
        //ftree[(convert_pop[*it_C3]*(convert_pop[*it_C3]-1.0))/2.0 + convert_pop[*it_C2]].mean_bc2a_c1++;
      }

    }
  }

}

void
update_jackknife(std::vector<threepop_stats>& ftree, double& mean_all, double& sd_all, int num_trees){

  double mean_block = 0.0;
  int k = 0;
  
  for(std::vector<threepop_stats>::iterator it_ftree = ftree.begin(); it_ftree != ftree.end(); it_ftree++){
    (*it_ftree).count_ac1_bc2 /= (double) num_trees;
    (*it_ftree).count_ac2_bc1 /= (double) num_trees;
    (*it_ftree).count_ac1b_c2 /= (double) num_trees;
    (*it_ftree).count_ac2b_c1 /= (double) num_trees;
    (*it_ftree).count_bc1a_c2 /= (double) num_trees;
    (*it_ftree).count_bc2a_c1 /= (double) num_trees;

    //(*it_ftree).mean_ac1_bc2  += (*it_ftree).count_ac1_bc2;
    //(*it_ftree).mean_ac2_bc1  += (*it_ftree).count_ac2_bc1;
    //(*it_ftree).mean_ac1b_c2  += (*it_ftree).count_ac1b_c2;
    //(*it_ftree).mean_ac2b_c1  += (*it_ftree).count_ac2b_c1;
    //(*it_ftree).mean_bc1a_c2  += (*it_ftree).count_bc1a_c2;
    //(*it_ftree).mean_bc2a_c1  += (*it_ftree).count_bc2a_c1;

    //(*it_ftree).sd_ac1_bc2  += (*it_ftree).count_ac1_bc2*(*it_ftree).count_ac1_bc2;
    //(*it_ftree).sd_ac2_bc1  += (*it_ftree).count_ac2_bc1*(*it_ftree).count_ac2_bc1;
    //(*it_ftree).sd_ac1b_c2  += (*it_ftree).count_ac1b_c2*(*it_ftree).count_ac1b_c2;
    //(*it_ftree).sd_ac2b_c1  += (*it_ftree).count_ac2b_c1*(*it_ftree).count_ac2b_c1;
    //(*it_ftree).sd_bc1a_c2  += (*it_ftree).count_bc1a_c2*(*it_ftree).count_bc1a_c2;
    //(*it_ftree).sd_bc2a_c1  += (*it_ftree).count_bc2a_c1*(*it_ftree).count_bc2a_c1;

    double tmp;
    tmp = (*it_ftree).count_ac1_bc2 - (*it_ftree).count_ac1b_c2 - (*it_ftree).count_bc1a_c2;
    (*it_ftree).mean_ac1bc2 += tmp;
    (*it_ftree).sd_ac1bc2   += tmp*tmp;
    mean_block              += tmp;

    tmp = (*it_ftree).count_ac2_bc1 - (*it_ftree).count_ac2b_c1 - (*it_ftree).count_bc2a_c1;
    (*it_ftree).mean_ac2bc1 += tmp;
    (*it_ftree).sd_ac2bc1   += tmp*tmp;
    mean_block              += tmp;

    (*it_ftree).count_ac1_bc2 = 0.0;
    (*it_ftree).count_ac2_bc1 = 0.0;
    (*it_ftree).count_ac1b_c2 = 0.0;
    (*it_ftree).count_ac2b_c1 = 0.0;
    (*it_ftree).count_bc1a_c2 = 0.0;
    (*it_ftree).count_bc2a_c1 = 0.0;
    k++;
  }

  mean_block /= k;
  mean_all += mean_block;
  sd_all   += mean_block*mean_block;

}

void
threepop(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("poplabels") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, poplabels, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Calculating 3 population ftree statistics for " << options["anc"].as<std::string>() << " .." << std::endl;

  Sample sample;
  sample.Read(options["poplabels"].as<std::string>());
  if(sample.groups.size() != 3){
    std::cerr << "poplabels file contains more than three populations" << std::endl;
    exit(1);
  }
  int pop1 = 0, pop2 = 1, pop3 = 2;

  int jackknife_bin = 500, num_blocks = 0, num_trees = 0;
  int a = 0, b = 0;

  std::vector<int> convert_pop(sample.group_of_haplotype.size(), -1), pop_a(sample.group_sizes[pop1], 0), pop_b(sample.group_sizes[pop2], 0), pop_c(sample.group_sizes[pop3], 0);
  int i1 = 0, i2 = 0, i3 = 0, k = 0;
  std::vector<int>::iterator it_samples = sample.group_of_haplotype.begin(), it_pop_a = pop_a.begin(), it_pop_b = pop_b.begin(), it_pop_c = pop_c.begin();
  for(std::vector<int>::iterator it_convert_pop = convert_pop.begin(); it_convert_pop != convert_pop.end(); it_convert_pop++){
    if(*it_samples == pop1){
      *it_convert_pop = i1;
      *it_pop_a = k;
      it_pop_a++;
      i1++;
    }
    if(*it_samples == pop2){
      *it_convert_pop = i2;
      *it_pop_b = k;
      it_pop_b++;
      i2++;
    }
    if(*it_samples == pop3){
      *it_convert_pop = i3;
      *it_pop_c = k;
      it_pop_c++;
      i3++;
    }
    k++;
    it_samples++;
  }

  ////////////
  
  std::vector<threepop_stats> ftree(sample.group_sizes[pop1]*sample.group_sizes[pop2]*(sample.group_sizes[pop3]*(sample.group_sizes[pop3]-1.0))/2.0);
  for(std::vector<threepop_stats>::iterator it_ftree = ftree.begin(); it_ftree != ftree.end(); it_ftree++){
    (*it_ftree).pop1 = sample.groups[pop1];
    (*it_ftree).pop2 = sample.groups[pop2];
    (*it_ftree).pop3 = sample.groups[pop3];
  }

  std::vector<threepop_stats>::iterator it_ftree = ftree.begin();
  for(int a1 = 0; a1 < pop_a.size(); a1++){
    for(int b1 = 0; b1 < pop_b.size(); b1++){
      for(int i1 = 0; i1 < pop_c.size(); i1++){
        for(int i2 = 0; i2 < i1; i2++){
          (*it_ftree).a  = pop_a[a1];
          (*it_ftree).b  = pop_b[b1];
          (*it_ftree).c1 = pop_c[i1];
          (*it_ftree).c2 = pop_c[i2];
          it_ftree++; 
        }
      }
    }
  }
  
  double mean_all = 0.0, sd_all = 0.0;  

  ///////

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
  std::vector<float> coords;
  float num_bases_tree_persists = 0.0;

  std::vector<int> tree_length(ancmut.NumTrees()), num_muts(ancmut.NumTrees());
  int tr = mut.info[0].tree;
  for(Muts::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
    if(tr == (*it_mut).tree){
      tree_length[tr] += (*it_mut).dist;
    }else{
      tr = (*it_mut).tree;
    }
  }
  tr = 0;
  while(num_bases_tree_persists >= 0.0){
    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    if(num_bases_tree_persists >= 0.0){
      num_muts[tr] = 0.0;
      for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
        num_muts[tr] += (*it_node).num_events;
      }
    }
    tr++;
  }
  ancmut.OpenFiles(options["anc"].as<std::string>(), options["mut"].as<std::string>());

  std::vector<int> tmp = tree_length;
  std::sort(tmp.begin(), tmp.end());
  int threshold_L = tmp[(int)(0.25*tmp.size())];
  tmp = num_muts;
  std::sort(tmp.begin(), tmp.end());
  int threshold_N = tmp[(int)(0.25*tmp.size())];

  //////////////////////////////////////
  int root = 2 * ancmut.NumTips() - 2;
  int N = ancmut.NumTips();
  num_bases_tree_persists = 0.0;
  num_trees = 0;
  tr = 0;
  while(num_bases_tree_persists >= 0.0){
    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

    if(num_bases_tree_persists >= 0.0){

      if(tree_length[tr] > threshold_L && num_muts[tr] > threshold_N){

        //do something with mtr
        mtr.tree.GetCoordinates(coords);

        //in principle, can iterate over a and b here
        //(and can make more efficient by )

        for(int a1 = 0; a1 < pop_a.size(); a1++){
          a = pop_a[a1];
          for(int b1 = 0; b1 < pop_b.size(); b1++){
            b = pop_b[b1];

            int mrca = -1, a_anc = a, b_anc = b, a_anc_prev = a, b_anc_prev = b;
            std::vector<int> ancestors_a(N+1), ancestors_b(N+1);
            int k1 = 1, k2 = 1;
            ancestors_a[0] = a;
            while(mtr.tree.nodes[a_anc].parent != NULL){
              a_anc = (*mtr.tree.nodes[a_anc].parent).label;
              ancestors_a[k1] = a_anc;
              k1++;
            }
            ancestors_b[0] = b;
            while(mtr.tree.nodes[b_anc].parent != NULL){
              b_anc = (*mtr.tree.nodes[b_anc].parent).label;
              ancestors_b[k2] = b_anc;
              k2++;
            }
            for(int i = 0; i < k1; i++){
              for(int j = 0; j < k2; j++){
                if(ancestors_a[i] == ancestors_b[j]){
                  mrca       = ancestors_a[i];
                  a_anc      = ancestors_a[i];
                  a_anc_prev = ancestors_a[i-1];
                  b_anc      = ancestors_b[j];
                  b_anc_prev = ancestors_b[j-1];
                  break;
                }
              }
              if(mrca != -1) break;
            } 
            assert(mrca != -1);

            std::vector<int> C1, C2, C3;

            get_C(mtr.tree.nodes[a_anc_prev], sample, pop3, C1);
            get_C(mtr.tree.nodes[b_anc_prev], sample, pop3, C2);

            //get C3 by looking at all Cs not in C1 or C2
            std::sort(C1.begin(), C1.end());
            std::sort(C2.begin(), C2.end());
            C3.resize(sample.group_sizes[pop3] - C1.size() - C2.size());

            if(C3.size() > 0){

              std::vector<int>::iterator it_C1 = C1.begin();
              std::vector<int>::iterator it_C2 = C2.begin();
              std::vector<int>::iterator it_C3 = C3.begin();
       
              int k = 0;
              for(std::vector<int>::iterator it_samples = sample.group_of_haplotype.begin(); it_samples != sample.group_of_haplotype.end(); it_samples++){
                if(*it_samples == pop3){ 
      
                  //std::cerr << k << " " << *it_C1 << " " << *it_C2 << std::endl; 
                  if(it_C1 != C1.end() && it_C2 != C2.end()){
                    if(*it_C1 == k){
                      it_C1++;
                    }else if(*it_C2 == k){
                      it_C2++;
                    }else{
                      *it_C3 = k;
                      it_C3++;              
                    }
                  }else if(it_C1 != C1.end()){            
                    if(*it_C1 == k){
                      it_C1++;
                    }else{
                      *it_C3 = k;
                      it_C3++;              
                    }
                  }else if(it_C2 != C2.end()){            
                    if(*it_C2 == k){
                      it_C2++;
                    }else{
                      *it_C3 = k;
                      it_C3++;              
                    } 
                  }else{
                    *it_C3 = k;
                    it_C3++;  
                  }
               
                }
                k++;
              }
            }

            //count topologies
            assign_ftree_3pop(a, b, C1, C2, C3, ftree, convert_pop, pop_a.size(), pop_b.size(), pop_c.size());

          }
        }

        //do jackknife

        //jackknife_bin
        //num_blocks
        num_trees++;

        if(num_trees == jackknife_bin){
          update_jackknife(ftree, mean_all, sd_all, num_trees);
          num_trees = 0; 
          num_blocks++;
        }

      }

      tr++;
    }

  }

  if(num_trees > 0){
    update_jackknife(ftree, mean_all, sd_all, num_trees);
    num_trees = 0; 
    num_blocks++;
  }


  mean_all /= num_blocks;
  sd_all    = sqrt( sd_all/(num_blocks * (num_blocks - 1.0)) - mean_all * mean_all/(num_blocks - 1.0) );

  for(std::vector<threepop_stats>::iterator it_ftree = ftree.begin(); it_ftree != ftree.end(); it_ftree++){

    (*it_ftree).mean_ac1bc2  /= num_blocks;
    (*it_ftree).mean_ac2bc1  /= num_blocks;
    
    (*it_ftree).mean_ac1_bc2 /= num_blocks;
    (*it_ftree).mean_ac2_bc1 /= num_blocks;
    (*it_ftree).mean_ac1b_c2 /= num_blocks;
    (*it_ftree).mean_ac2b_c1 /= num_blocks;
    (*it_ftree).mean_bc1a_c2 /= num_blocks;
    (*it_ftree).mean_bc2a_c1 /= num_blocks;

    (*it_ftree).sd_ac1bc2  = sqrt( (*it_ftree).sd_ac1bc2/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_ac1bc2 * (*it_ftree).mean_ac1bc2/(num_blocks - 1.0) ); 
    (*it_ftree).sd_ac2bc1  = sqrt( (*it_ftree).sd_ac2bc1/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_ac2bc1 * (*it_ftree).mean_ac2bc1/(num_blocks - 1.0) ); 

    (*it_ftree).sd_ac1_bc2 = sqrt( (*it_ftree).sd_ac1_bc2/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_ac1_bc2 * (*it_ftree).mean_ac1_bc2/(num_blocks - 1.0) );
    (*it_ftree).sd_ac2_bc1 = sqrt(  (*it_ftree).sd_ac2_bc1/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_ac2_bc1 * (*it_ftree).mean_ac2_bc1/(num_blocks - 1.0) ); 
    (*it_ftree).sd_ac1b_c2 = sqrt(  (*it_ftree).sd_ac1b_c2/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_ac1b_c2 * (*it_ftree).mean_ac1b_c2/(num_blocks - 1.0) );
    (*it_ftree).sd_ac2b_c1 = sqrt(  (*it_ftree).sd_ac2b_c1/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_ac2b_c1 * (*it_ftree).mean_ac2b_c1/(num_blocks - 1.0) );
    (*it_ftree).sd_bc1a_c2 = sqrt(  (*it_ftree).sd_bc1a_c2/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_bc1a_c2 * (*it_ftree).mean_bc1a_c2/(num_blocks - 1.0) );
    (*it_ftree).sd_bc2a_c1 = sqrt(  (*it_ftree).sd_bc2a_c1/(num_blocks*(num_blocks - 1.0)) - (*it_ftree).mean_bc2a_c1 * (*it_ftree).mean_bc2a_c1/(num_blocks - 1.0) );

  }

  std::ofstream os(options["output"].as<std::string>() + ".ftree");
  os << "pop1 pop2 pop3 mean_tf3 sd_tf3 num_blocks\n";
  os << pop1 << " " << pop2 << " " << pop3 << " " << mean_all << " " << sd_all << " " << num_blocks << "\n";
  os.close();

  if(0){
  std::ofstream os(options["output"].as<std::string>() + ".ftree");
  os << "pop1 pop2 pop3 a b c1 c2 mean_tf3 sd_tf3\n";
  for(std::vector<threepop_stats>::iterator it_ftree = ftree.begin(); it_ftree != ftree.end(); it_ftree++){
   
    os << (*it_ftree).pop1         << " " << (*it_ftree).pop2         << " " << (*it_ftree).pop3 << " ";
    os << (*it_ftree).a            << " " << (*it_ftree).b            << " " << (*it_ftree).c1 << " " << (*it_ftree).c2 << " ";
    os << (*it_ftree).mean_ac1bc2 << " " << (*it_ftree).sd_ac1bc2 << " ";
    os << "\n";

    os << (*it_ftree).pop1         << " " << (*it_ftree).pop2         << " " << (*it_ftree).pop3 << " ";
    os << (*it_ftree).a            << " " << (*it_ftree).b            << " " << (*it_ftree).c2 << " " << (*it_ftree).c1 << " ";
    os << (*it_ftree).mean_ac2bc1 << " " << (*it_ftree).sd_ac2bc1 << " ";
    os << "\n";

  } 
  os.close();
  }

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

