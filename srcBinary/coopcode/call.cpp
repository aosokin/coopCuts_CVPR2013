
#include <vector>
#include <cstdio>
#include "itBM.h"
#include "edat_read.h"
#include "costfun.h"

#include "time.h"

/*
  This example shows how to use the cooperative cut code. For the theory and 
  algorithm, see the paper

  S. Jegelka and J. Bilmes. "Submodularity beyond submodular energies: coupling
  edges in graph cuts". IEEE Conference on Computer Vision and Pattern Recognition 
  (CVPR), 2011.

  If you do not want to use OpenCV, you can replace the I/O yourself to use your 
  preferred functions or library.

  Anton Osokin: cut off OpenCV, make high order function piecewise-linear with one break point


  Edge classes are assumed to be in a different binary file that uses the 
  following format:
  number_of_classes
  class of edge 1
  class of edge 2 .....
  ... class of edge m  

  - 8-neighbor graph (incl. diagonals)
  - first set of edges: (vertical up)
    for (j=0; j<width; j++) {
      for (i=0; i<(height-1); i++) { 
      data[i*step + j*channels to data[(i+1)*step + j*channels]
  - second set: (horizontal to left)
    for (j=0; j<(width-1); j++) {
      for (i=0; i<(height); i++) {
       from data[i*step + j*channels + k]  to data[i*step + (j+1)*channels + k]
  - third set:
    for (j=0; j<(width-1); j++) {
      for (i=0; i<(height-1); i++) {
        data[i*step + j*channels + k] to data[(i+1)*step + (j+1)*channels + k]
  - fourth set:
    for (j=0; j<(width-1); j++) {
      for (i=0; i<(height-1); i++) {
        data[(i+1)*step + j*channels + k] to data[i*step + (j+1)*channels + k]


   The cut cost is 
   cost(C)= w(C \intersect E_terminal) + gamma * f(C intersect E_interpixel)
   where f(C) = \sum_i f_i(C \intersect S_i)
   and f_i(A) = w(A \intersect S_i)  if w(A \intersect S_i) <= theta*w(S_i);
    f_i(A) = theta*w(S_i) + sqrt(1+ w(A \intersect S_i) - theta*w(S_i)) -1 
    otherwise. Edges between nodes of identical color are treated in a separate
    group without any discounts.

  To use something else than sqrt, replacing the function 'sqrt2' in costfun.cpp
  should work (but no guarantees).
 */


int main (int argc, char** argv){

  if (argc < 6) {
    printf("usage: coopcode edgeClassFile edgeWeightsFile unaryFile lambda theta [outputLabelFile] [alpha]\n");
    return 1;
  }

  // parameters:
  double gamma = atof(argv[4]);
  double thresh = atof(argv[5]);
  char* classfile = argv[1];  // file that has the edge classes in binary format
  char* unaryFile = argv[3];  // file with the terminal weights
  char* edgeWeightsFile = argv[2]; // file with edge weights
  
  char* outputLabelFile = (argc >= 6)? argv[6] : NULL; // check if there is an output label file
  double alpha = (argc >= 7)? atof(argv[7]) : 0; // check if there exist an alpha

  /*
    1. create an edge data structure with edge weights and edge types/classes.
       This must be an Edge_data class and hold the weights, classes and image 
       size; Edat_read is one way to implement reading into this structure.
       You can also write your own, make it inherit from Edge_data.
  */
  Edat_read* edat = new Edat_read();
  int maxiter = 15;
  /*
    role of the modthresh parameter: In the CVPR paper's experiments, we kept
    one special class of edges that do not enjoy any discount. The edges in 
    that class are between very similar (distance between node color vectors smaller
    than modthresh).
   */
  double modthresh = 0.0000000001;

  /*
    read the weights form the weight file
  */
  if ( !edat -> readWeightsFromFile(edgeWeightsFile)) {
	  return 1;
  }
 

  /*
    read in the classes, the order of the classes must follow the order of the
    edges outlined above. For each edge described above, the following edge
    is the reverse edge, meaning if edge 2*i is (u,v), then edge 2*i+1 is (v,u).
   */
  if (! edat -> readClasses(classfile) ) {
	  return 1;
  }
  // alternatively, set the classes via an integer array of length 2*edat->m:
  // edat->readClasses(int* classes_in, int num_classes, double modthresh);
  // num_classes is the number of classes, modthresh determines
  // when two nodes are similar enough so that the edge between them is heavy
  // enough not to get discounts.

  /*
    2. Create a cost function and set gamma for it, and the unary weights.
       The unary waights must be an array of length 2*n, first the source, and then
       the sink weights. The order of the nodes is equivalent to the order
       in which opencv reads. To iterate through the nodes, use e.eg.
       for (int y = 0; y < edat->height; ++y) {
         for (int x = 0; x < edat->width; ++x) {
           nodeid = y*edat->width + x;
         }
       }
  */
  Costfun* Cf = new Costfun();
  Cf->setGamma(gamma);
  int n = edat->height * edat->width;

  Cf->setAlpha(alpha);

  /*
    read in unary weights
  */
  FILE* fp = fopen(unaryFile, "rb");
  if (fp == NULL) {
	  printf("Error while reading unary file: could not open file %s\n", unaryFile);
	  return 1;
  }
  double* unaryWeights = (double*) malloc(2 * n * sizeof(double));
  if ( fread(unaryWeights, sizeof(double), 2 * n, fp) != 2 * n) {
	  printf("Error while reading unary file: could not read 'unaryWeights': double[2 * numNodes]\n");
	  return 1;
  }

  Cf->setEdges(edat, unaryWeights, thresh);

  double tStart = (double)clock();

  /*
    3. Create the optimization algorithm. This requires a cost function,
       and the maximum number of iterations.
   */
  ItBM* algo = new ItBM(Cf, maxiter);

  std::vector<unsigned int> cut;
  bool sort_like_input = false;
  bool save_labels = false;

  /*
    4. Call the minimizer. "Cut" will in the end contain the set of cut edges,
       if save_labels is set, then we write out the node labels (as char)
       in binary format into the file labelfile (last argument).
  */
  double energy = algo->minimize(cut, sort_like_input, save_labels);

  printf("cut has %d edges\n", cut.size());

  // get the node labels
  std::vector<bool> nodelabels = algo->get_node_labels();

  double elapsedTime = (clock() - tStart) / CLOCKS_PER_SEC;
  printf("Time: %fs\n", elapsedTime);


  if (outputLabelFile != NULL) {
	 // write out node labels
	FILE* of = fopen(outputLabelFile, "wb");

	if (of == NULL) {
		printf("Could not create the output file\n.");
		return 1;

	}

	fwrite(&n, sizeof(int), 1, of);
	for(int iNode = 0; iNode < n; ++iNode){
		signed char curLabel = (nodelabels[iNode])? 1 : 0;
		fwrite(&curLabel, sizeof(signed char), 1, of);
	}

	// write elapsedTime to the file
	fwrite(&elapsedTime, sizeof(double), 1, of);

	// write the energy computed by the method
	fwrite(&energy, sizeof(double), 1, of);


	fclose(of);
  }



  /*
    show the image
   */  
  /*
  IplImage* img = cvLoadImage(imname);
  if(!img) { printf("Could not load image file!\n"); }
  uchar* data = (uchar *)img->imageData;
  int channels = img->nChannels, step = img->widthStep;
  int nodeid;

  for (int y = 0; y < edat->height; ++y) {
    for (int x = 0; x < edat->width; ++x) {
      nodeid = y*edat->width + x;
      if (!nodelabels[nodeid]) {
	for (int k=0; k<3; ++k)
	  data[y*step + x*channels + k] = 255;
      }
    }

  }

  cvNamedWindow("mainWin",CV_WINDOW_AUTOSIZE);
  cvMoveWindow("mainWin",100,100);
  cvShowImage("mainWin",img);
  cvWaitKey(0); 
  cvReleaseImage(&img);
  */
  delete algo;
  delete Cf;
  

  delete edat;

  return 0;  
}
