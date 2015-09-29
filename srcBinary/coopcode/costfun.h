#ifndef __COSTFUN_H__
#define __COSTFUN_H__

#include <vector>
#include <list>
#include "../graphCuts/maxflow-v3.02.src/graph.h"
#include "edge_data.h"

/*
  This class holds the graph, updates edge weights and computes cut costs.
  If you use this software for a publication, then cite the following papers:
  
  S. Jegelka and J. Bilmes. "Submodularity beyond submodular energies: coupling
  edges in graph cuts". IEEE Conference on Computer Vision and Pattern Recognition 
  (CVPR), 2011.
 */



class Costfun
{
 
 public:
  /*constructor and destructor:
   edge list defined by heads/tails should not contain back edge -- 
   the constructor will create the back edge (with the same weight)
   automatically. The classes array should be twice as long as weights, 
   where the classes of an edge and its sister come one after the other.
  */
  Costfun(double* weights, double* srcweights, double* sinkweights, int* classes, int* tails, int* heads, int m_, int nclasses_, int n_);
  Costfun();
  ~Costfun();

  Graph<double, double, double>* getGraph() { return myGraph; };

  // the algorithm will iteratively run maxflow and use S as new basis set
  double run_maxflow();
  void find_cut_set();

  void print_node_labels(Graph<double,double,double>::termtype default_label);

  /* computes sum for each class and whether it is satisfied,
     then updates edge weights accordingly to rho_e(E-e) or rho_e(S) 
  */
  // if edge indices refer to the internal sorted list elist (sorted by class)
  void reset_wts2(const std::vector<unsigned int>& S);
  // if edge indices refer to original input (order in heads/tails) 
  void reset_wts2_unsorted(std::vector<unsigned int>& S); // calls the above but changes S
  void reset_wtsCS();

  double subm_cut_cost(); // submodular cost of the current cut current_S
  void get_cut(std::vector<unsigned int>& cut);
  void get_original_indices(std::vector<unsigned int>& cut);
  std::vector<bool> get_node_labels(Graph<double,double,double>::termtype default_label);
  void get_node_labels(Graph<double,double,double>::termtype default_label, char* nodelabs);
  void setGamma(double gammaIn);

  int get_n() { return n; };

  void readFromImage(const char* imname, const char* unaryname,
		     double threshlim, double modthresh=1.0/1000.0, int nclass=10);
  void setEdges(Edge_data* Estuff, const double* unaryWeights, double threshlim);

  //AOSOKIN
  double getBaseFlow() const {return baseflow;}  

 private:
  unsigned int n, m;
  int nclasses;
  struct edgeInf;
  Graph<double, double, double>* myGraph;
  std::vector<edgeInf> elist;
  std::vector<unsigned int> csizes; //class sizes
  int lim;
  unsigned int submodlim; // pointer to first modular edge (class > lim)
  double gamma; // will influence terminal weights

  std::vector<double> rhoEe;
  std::vector<double> threshs;
  std::vector<unsigned int> current_S;
  // need to store source/sink weights to re-set capacities
  std::vector<double> termweights; // source - sink weight

  std::vector<unsigned int> sortindex;
  bool s_uptodate, cost_uptodate;
  double scost;

  /* 
     for each class, we pre-compute the weight of the class 
     and the cutoff-weight
  */
  std::vector<double> csums;
  std::vector<double> csumsTrunc;
  double baseflow;

  struct edgeInf
  {
    Graph<double,double,double>::arc_id edgep;
    unsigned int eindex;
    double weight;
    int classno;
  };

  friend bool comp(const Costfun::edgeInf&, const Costfun::edgeInf&);

  // AOSOKIN
  double sqrt2(double x);
  double alpha_;
public:
	void setAlpha(double alpha) { alpha_ = alpha;}


};



#endif
