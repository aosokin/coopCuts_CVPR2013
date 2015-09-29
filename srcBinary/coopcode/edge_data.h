#ifndef __EDGEDATA_H__
#define __EDGEDATA_H__

#include <cstdlib>

/*
  generic edge information class

  If you use this software for a publication, then cite the following papers:
  
  S. Jegelka and J. Bilmes. "Submodularity beyond submodular energies: coupling
  edges in graph cuts". IEEE Conference on Computer Vision and Pattern Recognition 
  (CVPR), 2011.
 */


class Edge_data {

 public:
  Edge_data() { 
    m=0; width=0; height=0; channels = 0;
    classes = NULL;
    weights = NULL;
  };

  ~Edge_data() {
    if (weights != NULL)
      delete[] weights;
    if (classes != NULL)
      delete[] classes;
  }

  double* weights;  // array of edge weights
  int* classes;  // array of edge classes
  unsigned int m, nclass;
  int height, width, channels;
  // m is the number of edges
  // nclass is the number of classes
  // channels is the number of channels

  // weights has length m, and classes has length 2*m, since edges
  // (i,j) and (j,i) have the same weights, but can have different
  // classes

};

#endif
