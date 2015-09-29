#ifndef __ITBM_H__
#define __ITBM_H__

#include <vector>
#include "costfun.h"

/*
  This part implements the iterative minimization, and takes any costfunction object.
  If you use this software for a publication, then cite the following papers:
  
  S. Jegelka and J. Bilmes. "Submodularity beyond submodular energies: coupling
  edges in graph cuts". IEEE Conference on Computer Vision and Pattern Recognition 
  (CVPR), 2011.
 */


class ItBM
{

 public:

  ItBM(Costfun* Fin, int maxiter_);

  // minimization:
  // second bool tells whether to save node labels
  double minimize(std::vector<unsigned int>& cut, bool, bool, const char* outfile=NULL);

  std::vector<bool> get_node_labels();
  
  void append_blabels(const char* outfile); // writes number of nodes and labels into a file, appends if the file exists
  void append_blabels(const char* outfile, double* fval); 

 private:

  Costfun* F;
  int maxiter;
  bool is_sorted;
  std::vector<bool> best_node_labels;

  void append_thoselabels(const char* outfile, std::vector<bool> thoselabels, double* fval);

};

#endif
