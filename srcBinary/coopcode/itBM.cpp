#include <iostream>
#include <stdexcept>
#include <sys/stat.h>
#include <cstdio>
#include "itBM.h"
#include "costfun.h"


ItBM::ItBM(Costfun* Fin, int maxiter_) {

  F = Fin;
  maxiter = maxiter_;
  is_sorted = false;

}

/*
  minimize. If save_labels=true, then the foudn node labeld are
  saved in the file indicated by outfile
  The numbers of the cut edges will be in 'cut'.
  If sort_like_input is set, then the edge numbers refer to the input order,
  otherwise to the re-sorted order (resorting was done by class labels)
*/
double ItBM::minimize(std::vector<unsigned int>& cut, bool sort_like_input, bool save_labels, const char* outfile) {

  if (save_labels && outfile==NULL)
    save_labels = false;

  // first iteration: empty set
  double currentmin_submod, newobj, min2;
  int iter = 0;

  double tmp  = F -> run_maxflow();

  F -> find_cut_set();
  newobj = F->subm_cut_cost();
  currentmin_submod = newobj + 1; // just ensure that it's bigger

  while ((newobj < currentmin_submod) && (iter < maxiter)) {

    iter++;

    currentmin_submod = newobj;
    printf("current min: %f ", currentmin_submod);
    F -> get_cut(cut);
    best_node_labels = F->get_node_labels(Graph<double,double,double>::SOURCE);

    // re-iterate
    F -> reset_wtsCS();
    F -> run_maxflow();
    F -> find_cut_set();
    newobj = F->subm_cut_cost();
    printf(" newobj=%f\n", newobj);    
  }
  
  if (save_labels) {
    append_blabels(outfile, &currentmin_submod);
  }
  
  if (sort_like_input) {
    F -> get_original_indices(cut);
  }
  return currentmin_submod + F -> getBaseFlow();

}




std::vector<bool> ItBM::get_node_labels() {

  if (best_node_labels.size()>0)
    return best_node_labels;
  else
    return F->get_node_labels(Graph<double,double,double>::SOURCE);

}


/* 
   check whether a file exists
   taken from http://www.techbytes.ca/techbyte103.html
 */

bool fileExists(const char* strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename,&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}


void ItBM::append_blabels(const char* outfile){

  int n = best_node_labels.size();
  char* nodelabs = (char*) malloc (n * sizeof(char));
  for (int i=0; i<n; i++){
    nodelabs[i] = best_node_labels[i];
  }
  FILE* fp;
  if (fileExists(outfile)) {
    fp = fopen(outfile,"a+b");
  } else {
    printf("starting file\n");
    fp = fopen(outfile,"wb");
    fwrite(&n, sizeof(int), 1, fp);
  }
  
  fwrite(nodelabs, sizeof(char), n, fp);
  fclose(fp);
  
  free(nodelabs);

}


void ItBM::append_blabels(const char* outfile, double* fval){

  int n = best_node_labels.size();
  char* nodelabs = (char*) malloc (n * sizeof(char));
  for (int i=0; i<n; i++){
    nodelabs[i] = best_node_labels[i];
  }
  FILE* fp;
  if (fileExists(outfile)) {
    fp = fopen(outfile,"a+b");
  } else {
    printf("starting file\n");
    fp = fopen(outfile,"wb");
    fwrite(&n, sizeof(int), 1, fp);
  }

  fwrite(fval, sizeof(double), 1, fp);
  fwrite(nodelabs, sizeof(char), n, fp);
  fclose(fp);
  
  free(nodelabs);

}



