#ifndef TREE_HPP
#define TREE_HPP

////////////////////////////
// Class for Trees and AncesTrees
///////////////////////////

#include "gzstream.hpp"
#include "data.hpp"
#include "sample.hpp"

#include <iostream>
#include <iomanip>
#include <list>
#include <limits>
#include <deque>
#include <ctime>
#include <tgmath.h>

//Data structure of a node in a tree
struct Node{

  Node* parent = NULL;
  Node* child_left = NULL;
  Node* child_right = NULL;
  int label;

  float num_events = 0.0; //on branch on top of parent
  int SNP_begin = 0;
  int SNP_end = 0;
  double branch_length = 0.0; //of branch on top of parent

  void operator=(const Node& node){
    parent        = node.parent;
    child_left    = node.child_left;
    child_right   = node.child_right;
    label         = node.label;
    num_events    = node.num_events;
    SNP_begin     = node.SNP_begin;
    SNP_end       = node.SNP_end;
    branch_length = node.branch_length;
  }
  bool operator==(Node node) const {
    if(node.label==label){
      return true;
    }else{
      return false;
    }
  }
  bool operator!=(Node node) const {
    if(node.label!=label){
      return true;
    }else{
      return false;
    }
  }

};

//Data structure for leaves below a branch
struct Leaves{

  int num_leaves;
  std::vector<int> member;

};

//Class for trees
class Tree{

  private:

    void TraverseTreeToGetCoordinates(Node& n, std::vector<float>& coordinates);
    void TraverseTreeToGetCoordinates_sample_age(Node& n, std::vector<float>& coordinates);
    void GetSubTree(std::vector<int>& subpop, Tree& subtree, std::vector<int>& convert_index, std::vector<int>& number_in_subpop) const;

  public:

    std::vector<double>* sample_ages = NULL;
    std::vector<Node> nodes;

    //construct tree
    Tree(){};
    Tree(const Tree& tr){
      nodes.resize(tr.nodes.size());
      for(int i = 0; i < (int) nodes.size(); i++){
        if(tr.nodes[i].parent != NULL) nodes[i].parent = &nodes[(*tr.nodes[i].parent).label];
        else nodes[i].parent = NULL;
        if(tr.nodes[i].child_left != NULL) nodes[i].child_left = &nodes[(*tr.nodes[i].child_left).label];
        else nodes[i].child_left = NULL;
        if(tr.nodes[i].child_right != NULL) nodes[i].child_right = &nodes[(*tr.nodes[i].child_right).label];
        else nodes[i].child_right = NULL;
        nodes[i].label = tr.nodes[i].label;
        nodes[i].num_events = tr.nodes[i].num_events;
        nodes[i].branch_length = tr.nodes[i].branch_length;
        nodes[i].SNP_begin     = tr.nodes[i].SNP_begin;
        nodes[i].SNP_end       = tr.nodes[i].SNP_end;
      }
    };
    Tree(const Tree& tr, std::vector<double>& i_sample_ages){
      nodes.resize(tr.nodes.size());
      for(int i = 0; i < (int) nodes.size(); i++){
        if(tr.nodes[i].parent != NULL) nodes[i].parent = &nodes[(*tr.nodes[i].parent).label];
        else nodes[i].parent = NULL;
        if(tr.nodes[i].child_left != NULL) nodes[i].child_left = &nodes[(*tr.nodes[i].child_left).label];
        else nodes[i].child_left = NULL;
        if(tr.nodes[i].child_right != NULL) nodes[i].child_right = &nodes[(*tr.nodes[i].child_right).label];
        else nodes[i].child_right = NULL;
        nodes[i].label = tr.nodes[i].label;
        nodes[i].num_events = tr.nodes[i].num_events;
        nodes[i].branch_length = tr.nodes[i].branch_length;
        nodes[i].SNP_begin     = tr.nodes[i].SNP_begin;
        nodes[i].SNP_end       = tr.nodes[i].SNP_end;
      }
      sample_ages = &i_sample_ages;
    };
 
    void ReadTree(const char* line, int N);
    void WriteNewick(const std::string& filename_newick, double factor, const bool add = 0) const; //this is slow. For fast output use WriteOrientedTree
    void WriteNHX(const std::string& filename_nhx, std::vector<std::string>& property, const bool add = 0) const; 

    void FindAllLeaves(std::vector<Leaves>& leaves) const;
    void FindLeaves(Node& node, std::vector<Leaves>& leaves) const; //recursive algorithm to find leaves. stored in leaves.

    void GetCoordinates(std::vector<float>& coordinates);
    void GetCoordinates(int node, std::vector<float>& coordinates);

    void GetNumberOfLeavesInSubpop(const Node& n, std::vector<int>& subpop, std::vector<int>& number_in_subpop) const; 
    void GetSubTree(Sample& sample, Tree& subtree) const;
    void GetSubTree(Sample& sample, Tree& subtree, std::vector<int>& convert_index, std::vector<int>& number_in_subpop) const;

    void operator=(const Tree& tr){
      nodes.resize(tr.nodes.size());
      for(int i = 0; i < (int) nodes.size(); i++){
        if(tr.nodes[i].parent != NULL) nodes[i].parent = &nodes[(*tr.nodes[i].parent).label];
        else nodes[i].parent = NULL;
        if(tr.nodes[i].child_left != NULL) nodes[i].child_left = &nodes[(*tr.nodes[i].child_left).label];
        else nodes[i].child_left = NULL;
        if(tr.nodes[i].child_right != NULL) nodes[i].child_right = &nodes[(*tr.nodes[i].child_right).label];
        else nodes[i].child_right = NULL;
        nodes[i].label = tr.nodes[i].label;
        nodes[i].num_events = tr.nodes[i].num_events;
        nodes[i].branch_length = tr.nodes[i].branch_length;
        nodes[i].SNP_begin     = tr.nodes[i].SNP_begin;
        nodes[i].SNP_end       = tr.nodes[i].SNP_end;
      }
    }

};

//Data Structure for a mancinal tree along the genome
struct MarginalTree{
  int pos;
  Tree tree;

  MarginalTree(){}
  MarginalTree(int pos, Tree tree): pos(pos), tree(tree){}
  
  void Read(const std::string& line, int N);
  void Read(const std::string& line, int N, std::vector<double>& sample_ages);
  void Dump(FILE *pfile);
  void operator=(const MarginalTree& mtr){
    pos = mtr.pos;
    tree = mtr.tree;
  }

};

typedef std::list<MarginalTree> CorrTrees; //I could change this to deque but deque has no splice

//Used to find equivalent nodes across trees
struct EquivalentNode{

  int node1;
  int node2;
  float corr;

  EquivalentNode(){};
  EquivalentNode(int node1, int node2, float corr): node1(node1), node2(node2), corr(corr){};

  bool operator > (const EquivalentNode& n) const{
    return (corr > n.corr);
  }

};

//Calculates the Pearson correlation between two sets
class Correlation{

  private:

    int N;
    float N_float;

  public:

    Correlation(int N): N(N){
      N_float = (float) N;
    }
    float Pearson(const Leaves& set1, const Leaves& set2);

};

//Data structure for AncesTrees/Tree sequences
class AncesTree{

  private:

  public:

    std::vector<double> sample_ages;
    CorrTrees seq; //starting position on the genome 
    int N, L;

    void Read(igzstream& is);
    void Read(const std::string& filename); //read anc in long-format. 

    void Dump(FILE* pfile);
    void Dump(std::ofstream& os); 
    void Dump(const std::string& filename); //dump anc in long-format

    //Read ancs from other filetypes
    void ReadArgweaverSMC(const std::string& filename);
    void ReadRent(const std::string& filename, float Ne);
    void ReadNewick(const std::string& filename, float Ne);

    //Associate equivalent branches in seq
    void BranchAssociation(const Tree& ref_tree, const Tree& tree, std::vector<int>& equivalent_branches, std::vector<std::vector<int>>& potential_branches, int N, int N_total, float threshold_brancheq);
    void AssociateEquivalentBranches();

};

#endif //TREE_HPP
