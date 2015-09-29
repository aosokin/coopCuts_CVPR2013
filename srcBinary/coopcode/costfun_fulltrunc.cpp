#include "costfun.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <numeric>
#include <cstdio>
#include "../graphCuts/maxflow-v3.02.src/graph.h"
#include "edge_data.h"
#include "edat_read.h"

// comparator function using class numbers

bool comp(const Costfun::edgeInf& a, const Costfun::edgeInf& b)
{
  return a.classno < b.classno;
}



Costfun::Costfun(): n(0), m(0), nclasses(0), myGraph(NULL),lim(0), submodlim(0), 
		    gamma(1.0), s_uptodate(false), cost_uptodate(false), scost(0),
		    baseflow(0) {
  
  
}



Costfun::~Costfun(){

  delete myGraph;

}



void Costfun::setGamma(double gammaIn) {

  gamma = gammaIn;

  if (myGraph != NULL && n>0) {

    // re-set terminal weights

    double tg = 1/gamma;
    for (unsigned int i=0; i<n; i++){

      myGraph -> set_trcap(i, tg*termweights[i]);
      
    }
    cost_uptodate = false;
    
  }
}





void Costfun::reset_wts2_unsorted(std::vector<unsigned int>& S) {

  for (unsigned int i=0; i<S.size(); i++){
    S[i] = sortindex[ S[i] ];
  }

  std::sort(S.begin(), S.end());
  reset_wts2(S);

}

// THIS HAS CHANGED FOR FULL TRUNCATION
double Costfun::sqrt2(double x) {
  // old: return (sqrt(x+1)-1);

	// AOSOKIN
	return alpha_ * x;
  //return 0;
}


/*
  reset weights with respect to reference set S
*/

void Costfun::reset_wts2(const std::vector<unsigned int>& S)
{

  // assume S is ordered as elist, and entries in S are indices in elist
  
  std::vector<double> clws(nclasses,0.0); // to store class weight sums

  unsigned int i;
  int cl;
  double tmp;
  unsigned int lastone = S.size();

  for (i=0; i<S.size(); i++){
    cl = elist[ S[i]].classno;

    if (cl >= lim) {
      lastone = i;
      break;
    }
      
    clws[ cl ] += elist[ S[i]].weight;
    
    // set the rhoe if that hasn't happened yet
    if (rhoEe[ S[i]] < 0){
      tmp = (csums[cl] - elist[S[i]].weight); // weight of the other edges
      if ( tmp < threshs[cl]) {
        rhoEe[ S[i]] = csumsTrunc[cl] - tmp;
      } else {
        rhoEe[ S[i]] = csumsTrunc[cl] - ( threshs[cl] + sqrt2( tmp - threshs[cl]) );
      }
    }
  }
  
  // now we know class weights
  // so go through elist and re-set weights

  unsigned int counter = 0;

  for (unsigned int iC=0; iC<(csizes.size()-1); iC++){

    // is class saturated?
    if (clws[iC] >= threshs[iC]) {

      double base1 = sqrt2(clws[iC] - threshs[iC]);
      double base2 = clws[iC] - threshs[iC];

      for (i=0; i<csizes[iC]; i++){
        tmp = sqrt2( elist[ counter + i ].weight + base2) - base1;
        myGraph -> set_rcap(elist[i + counter].edgep, tmp);
      }

    } else {

      double miss = threshs[iC] - clws[iC];
      
      for (i=0; i<csizes[iC]; i++){
        
        if (elist[ counter + i].weight <= miss){
          tmp = elist[ counter + i].weight;
        } else {
          tmp = miss + sqrt2( elist[ counter + i].weight - miss );
        }
        myGraph -> set_rcap(elist[i + counter].edgep, tmp);
      }
    }
    counter += csizes[iC];    
  }
  
  for (i=0; i<lastone; i++){
    myGraph -> set_rcap( elist[ S[i]].edgep, rhoEe[ S[i]]);
  }

  //re-set modular and terminal edges
  for (i=submodlim; i<elist.size(); i++) {
    myGraph -> set_rcap( elist[i].edgep, elist[i].weight);
  }
  double tg = 1/gamma;
  for (i=0; i<n; i++){
    myGraph -> set_trcap(i, tg*termweights[i]);
  }
  
}



void Costfun::reset_wtsCS() {

  reset_wts2( current_S );

}


void Costfun::find_cut_set(){
  // identifies the currently cut edges
  // uses maxflow-given labels

  if (s_uptodate) { 
    return; 
  }

  //current_S = std::vector<unsigned int>(0);
  current_S.resize(0, 0);
  unsigned int i;
  double rcap = 0;
  Graph<double,double,double>::node_id head, tail;

  for (i=0; i < elist.size(); i++){
    rcap = myGraph -> get_rcap( elist[i].edgep );
    if (rcap == 0) {
      // get node capacities
	  myGraph -> get_arc_ends(elist[i].edgep, head, tail);
      if (myGraph -> what_segment(head,Graph<double,double,double>::SOURCE) == Graph<double,double,double>::SOURCE && myGraph -> what_segment(tail,Graph<double,double,double>::SOURCE) == Graph<double,double,double>::SINK) {
        current_S.push_back(i);
	  }
    }
  }
  s_uptodate = true;

}



double Costfun::run_maxflow() {

  s_uptodate = false;
  cost_uptodate = false;
  return myGraph -> maxflow();

} 




void Costfun::get_cut(std::vector<unsigned int>& cut) {

  if (!s_uptodate) {
    find_cut_set();
  }
  cut = current_S;

}


double Costfun::subm_cut_cost() {

  if (cost_uptodate)
    return scost;

  if (!s_uptodate) {
    find_cut_set();
  }

  std::vector<double> cls(nclasses+1,0);

  for (unsigned int i=0; i<current_S.size(); i++) {

    cls[ elist[ current_S[i]].classno ] += elist[ current_S[i]].weight;

  }

  double ret = 0;
  double tmpx = 0;
  for (int i=0; i<nclasses; i++) {
	  tmpx += cls[i];
    if (cls[i] < threshs[i]) {
      ret += cls[i];
    } else {
		double tmp = threshs[i] + sqrt2( cls[i] - threshs[i] );
      ret += threshs[i] + sqrt2( cls[i] - threshs[i] );
    }
  }
  ret += cls[nclasses];

  double tg = 1/gamma;
  for (unsigned int i=0; i<n; i++){
    if (myGraph -> what_segment(i,Graph<double,double,double>::SOURCE)==Graph<double,double,double>::SOURCE) {
      // node is in source segment: cut edge to sink
      
      ret += ( (termweights[i]>0) ? 0 : (-tg*termweights[i]) );

    } else {
      // node is in sink segment: cut cost max(termweights[i],0)
      ret += ( (termweights[i]>0) ? tg*termweights[i] : 0);
    }
  }

  scost = ret;
  cost_uptodate = true;

  return ret;

}


// get indices like in the input heads/tails, and not like the sorted
// edge list
void Costfun::get_original_indices(std::vector<unsigned int>& cut) {
  for (unsigned int i=0; i<cut.size(); i++){
    cut[i] = (elist[i]).eindex;
  }
}


std::vector<bool> Costfun::get_node_labels(Graph<double,double,double>::termtype default_label) {

  std::vector<bool> nodelabs(n);
  for (unsigned int i=0; i<n; i++) {
    nodelabs[i] = ( myGraph->what_segment(i,default_label) ==  Graph<double,double,double>::SOURCE);
  }

  return nodelabs;

}


void Costfun::get_node_labels(Graph<double,double,double>::termtype default_label, char* nodelabs) {

  for (unsigned int i=0; i<n; i++) {
    nodelabs[i] = ( myGraph->what_segment(i,default_label) ==  Graph<double,double,double>::SOURCE);
  }

}






// unaryWeights is array of length 2n: source weights, then sink weights

void Costfun::setEdges(Edge_data* Estuff, const double* unaryWeights,
		       double threshlim) {

  n = Estuff->height * Estuff->width;
  m = Estuff -> m;
  nclasses = Estuff->nclass;

  // add terminal edges
  myGraph = new Graph<double, double, double>(n+1,2*m+1);
  myGraph -> add_node(n);
  elist = std::vector<edgeInf>(2*m);

  sortindex = std::vector<unsigned int>(2*m);
  cost_uptodate = false;

  //put the weights into the termweights vector
  unsigned int i, j;
  baseflow = 0;
  double tg = 1/gamma;
  termweights = std::vector<double>(n);
  for (i=0; i<n; i++){
    
    termweights[i] = (unaryWeights[i] - unaryWeights[n+i]);
    if (termweights[i] > 0) {
      myGraph -> add_tweights(i, tg*termweights[i], 0);
      baseflow += tg * unaryWeights[n+i];   //BUG
    } else {
      myGraph -> add_tweights(i, 0, -tg*termweights[i]);
      baseflow += tg * unaryWeights[i]; //BUG
    }
  }
  
  printf("baseflow: %f\n", baseflow);

  // add interpixel edges
  csizes = std::vector<unsigned int>(nclasses+1,0);
  threshs = std::vector<double>(nclasses,0.0);
  csums = std::vector<double>(nclasses,0.0);
  csumsTrunc = std::vector<double>(nclasses);
  lim = nclasses;
  
  std::vector<double> cmaxs(nclasses,0);
  Graph<double,double,double>::arc_id farc = myGraph -> get_first_arc();

  // add edges to graph, row-major node numbering
  // (i,j)  is i*width+j
  // vertical
  unsigned int ct = 0;
  for (j=0; (int)j< Estuff->width; j++) {
    for (i=0; (int)i< (Estuff->height)-1; i++) {

      edgeInf nedge;
      myGraph -> add_edge( i*Estuff->width + j, (i+1)*Estuff->width + j, 50.0*Estuff->weights[ct], 50.0*Estuff->weights[ct]);
      nedge.edgep = farc++;
      nedge.eindex = 2*ct;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct];
      if (nedge.classno > nclasses) {
	printf("nclasses=%d, classno(%d)=%d\n", nclasses, 2*ct, nedge.classno);
      }

      if (nedge.classno >= lim){
	nedge.classno = lim;
      } else {
	csums[nedge.classno] += nedge.weight;
	if (nedge.weight > cmaxs[nedge.classno]){
	  cmaxs[nedge.classno] = nedge.weight;
	}
      }
      elist[2*ct] = nedge;
      csizes[nedge.classno]++;

      //sister edge
      nedge.edgep = farc++;
      nedge.eindex = 2*ct+1;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct + 1];
      if (nedge.classno >= lim){
        nedge.classno = lim;
      } else {
        csums[nedge.classno] += nedge.weight;
        if (nedge.weight > cmaxs[nedge.classno]){
          cmaxs[nedge.classno] = nedge.weight;
        }
      }
      elist[2*ct+1] = nedge;
      csizes[nedge.classno]++;

      ct++;
    }
  }

  // horizontal
  for (j=0; (int)j< (Estuff->width)-1; j++) {
    for (i=0; (int)i< (Estuff->height); i++) {

      edgeInf nedge;
      myGraph -> add_edge( i*Estuff->width + j+1, i*Estuff->width + j, 50.0*Estuff->weights[ct], 50.0*Estuff->weights[ct]);
      nedge.edgep = farc++;
      nedge.eindex = 2*ct;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct];
      if (nedge.classno >= lim){
	nedge.classno = lim;
      } else {
	csums[nedge.classno] += nedge.weight;
	if (nedge.weight > cmaxs[nedge.classno]){
	  cmaxs[nedge.classno] = nedge.weight;
	}
      }
      elist[2*ct] = nedge;
      csizes[nedge.classno]++;

      //sister edge
      nedge.edgep = farc++;
      nedge.eindex = 2*ct+1;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct + 1];
      if (nedge.classno >= lim){
        nedge.classno = lim;
      } else {
        csums[nedge.classno] += nedge.weight;
        if (nedge.weight > cmaxs[nedge.classno]){
          cmaxs[nedge.classno] = nedge.weight;
        }
      }
      elist[2*ct+1] = nedge;
      csizes[nedge.classno]++;

      ct++;
    }
  }


  // right diag down
  for (j=0; (int)j< (Estuff->width)-1; j++) {
    for (i=0; (int)i< (Estuff->height)-1; i++) {
      
      edgeInf nedge;
      myGraph -> add_edge( i*Estuff->width + j, (i+1)*Estuff->width + j+1, 50.0*Estuff->weights[ct], 50.0*Estuff->weights[ct]);
      nedge.edgep = farc++;
      nedge.eindex = 2*ct;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct];
      if (nedge.classno >= lim){
	nedge.classno = lim;
      } else {
	csums[nedge.classno] += nedge.weight;
	if (nedge.weight > cmaxs[nedge.classno]){
	  cmaxs[nedge.classno] = nedge.weight;
	}
      }
      elist[2*ct] = nedge;
      csizes[nedge.classno]++;

      //sister edge
      nedge.edgep = farc++;
      nedge.eindex = 2*ct+1;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct + 1];
      if (nedge.classno >= lim){
        nedge.classno = lim;
      } else {
        csums[nedge.classno] += nedge.weight;
        if (nedge.weight > cmaxs[nedge.classno]){
          cmaxs[nedge.classno] = nedge.weight;
        }
      }
      elist[2*ct+1] = nedge;
      csizes[nedge.classno]++;

      ct++;
    }
  }


  // left diag down
  for (j=0; (int)j< (Estuff->width)-1; j++) {
    for (i=0; (int)i< (Estuff->height)-1; i++) {

      edgeInf nedge;
      myGraph -> add_edge( (i+1)*Estuff->width + j, i*Estuff->width + j+1, 50.0*Estuff->weights[ct], 50.0*Estuff->weights[ct]);
      nedge.edgep = farc++;
      nedge.eindex = 2*ct;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct];
      if (nedge.classno >= lim){
	nedge.classno = lim;
      } else {
	csums[nedge.classno] += nedge.weight;
	if (nedge.weight > cmaxs[nedge.classno]){
	  cmaxs[nedge.classno] = nedge.weight;
	}
      }
      elist[2*ct] = nedge;
      csizes[nedge.classno]++;

      //sister edge
      nedge.edgep = farc++;
      nedge.eindex = 2*ct+1;
      nedge.weight = 50.0*Estuff->weights[ct];
      nedge.classno = Estuff->classes[2*ct + 1];
      if (nedge.classno >= lim){
        nedge.classno = lim;
      } else {
        csums[nedge.classno] += nedge.weight;
        if (nedge.weight > cmaxs[nedge.classno]){
          cmaxs[nedge.classno] = nedge.weight;
        }
      }
      elist[2*ct+1] = nedge;
      csizes[nedge.classno]++;

      ct++;
    }
  }

  submodlim = 2*m - csizes[csizes.size()-1];
  std::stable_sort(elist.begin(), elist.end(), ::comp);
  for (i=0; i<2*m; i++){
    sortindex[elist[i].eindex] = i;
  }
  
  double tmp;
  for (i=0; (int)i<nclasses; i++){
    tmp = threshlim*csums[i]; // threshold starts at threshlim fraction of class weight
    
    if (tmp >= cmaxs[i]){
      threshs[i] = tmp;
    } else {
      threshs[i] = cmaxs[i];
    }
    csumsTrunc[i] = threshs[i] + sqrt2(csums[i] - threshs[i]);
    
  }  

  rhoEe = std::vector<double>(2*m, -1);

}


/*
  average error, i.e. percentage of nodes where inlab and truelab disagree
*/
double comperr(char* inlab, char* truelab, int n) {

  double ret = 0.0;
  
  for (int i=0; i<n; i++) {

    if (inlab[i] != truelab[i])
      ret += 1;
  }

  ret /= (double)n;
  return ret;

}


