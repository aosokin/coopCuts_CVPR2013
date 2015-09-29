#include "edge_data.h"
#include "edat_read.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <cstdlib>

// destructor
Edat_read::~Edat_read() {
  free(weights);
  free(classes);
  weights = NULL;
  classes = NULL;
}


bool Edat_read::readWeightsFromFile(const char* edgeWeightsFile)
{
	bool resultFlag = true;
	if (edgeWeightsFile != NULL){
		try{
			// read the edge weights (created by the matlab code, makeweights.m)
			FILE* fp = fopen(edgeWeightsFile, "rb");
			if(fp == NULL) { 
				throw "could not open the file";
			}
			int height, width, weightsNumber;
			if (fread(&height, sizeof(int), 1, fp) != 1) { 
				fclose(fp);
				throw "could not read 'height': int";
			}
			if (fread(&width, sizeof(int), 1, fp) != 1) { 
				fclose(fp); 
				throw "could not read 'width': int"; 
			}
			if (fread(&weightsNumber, sizeof(int), 1, fp) != 1) { 
				fclose(fp); 
				throw "could not read 'weightsNumber': int"; 
			}
			double* weights = (double*) malloc(weightsNumber * sizeof(double));
			if (fread(weights, sizeof(double), weightsNumber, fp) != weightsNumber) {
				fclose(fp); 
				free(weights); 
				throw "could not read 'weights': double[weightsNumber]"; 
			}
			fclose(fp);
			fp = NULL;

			
			this -> height = height;
			this -> width = width;
			this -> m = weightsNumber;
			this -> channels = 3;

			if (weightsNumber != 4*height*width - 3*height - 3*width + 2) {
				free(weights);
				throw "the number of weights is incorrect: has to be 4*height*width - 3*height - 3*width + 2";
			}

			
			if (this -> weights != NULL) {
				free(this -> weights);
			}
			this -> weights = weights;
			weights = NULL;

			// printf("Successfully substituted the weight file\n");
		}
		catch( char * str ) {
			printf("Error while reading edgeWeightsFile: %s \n", str);
			resultFlag = false;
		}
	} else
	{
		return false;
	}
	return resultFlag;
}


/*
// The initial Stefanie's function
// create edge weights from an input image
void Edat_read::readWeightsFromImage(const char* imagename) {

  IplImage* img = cvLoadImage(imagename);
  if(!img) { printf("Could not load image file!\n"); }

  uchar *data;
  int step;
  int i,j,k;
  
  data      = (uchar *)img->imageData;
  height    = img->height;
  width     = img->width;
  step      = img->widthStep;
  channels  = img->nChannels;
  
  m = 4*height*width - 3*height - 3*width + 2;

  double tmp;
  unsigned int ei = 0;
  double sigma = 0.0;  

  weights = (double*) malloc(m*sizeof(double));
  
  // first set of edges
  for (j=0; j<width; j++) {
    for (i=0; i<(height-1); i++) {          
      weights[ei] = 0.0;      
      for (k=0; k<channels; k++) {
	tmp = double(data[(i+1)*step + j*channels + k]) - double(data[i*step + j*channels + k]) ;
	weights[ei] += tmp*tmp;
      }
      sigma += weights[ei]/m;
      ei++;
    }
  }
  
  // second set
  for (j=0; j<(width-1); j++) {
    for (i=0; i<(height); i++) {
      weights[ei] = 0.0;
      for (k=0; k<channels; k++) {	
	tmp = - (double)data[i*step + (j+1)*channels + k] + (double)data[i*step + j*channels + k];
	weights[ei] +=  tmp*tmp;
      }
      sigma += weights[ei]/m;
      ei++;
    }
  }
  
  // third set (diagonal)
  for (j=0; j<(width-1); j++) {
    for (i=0; i<(height-1); i++) {
      weights[ei] = 0.0;
      for (k=0; k<channels; k++) {
	tmp = (double)data[(i+1)*step + (j+1)*channels + k] - (double)data[i*step + j*channels + k];
	weights[ei] += tmp*tmp;
      }
      sigma += weights[ei]/m;
      ei++;
    }
  }

  // fourth set
  for (j=0; j<(width-1); j++) {
    for (i=0; i<(height-1); i++) {
      weights[ei] = 0.0;

      for (k=0; k<channels; k++) {
	tmp = - (double)data[(i+1)*step + j*channels + k] + (double)data[i*step + (j+1)*channels + k];
	weights[ei] +=  tmp*tmp;
      }
      sigma += weights[ei]/m;
      ei++;
    }
  }

  cvReleaseImage(&img);
  
  for (i=0; i<(int)m; i++) {
    weights[i] = 0.05 + 0.95*exp(-weights[i]/(2*sigma));
  }

  printf("sigma = %1.3e\n", sigma);

}

*/



// read in the edge classes from a binary file
// if you do not want modular treatment, then set modthresh < 0
// classfile should be of the form:
// num_classes
// int array of edge classes, length 2*m
bool Edat_read::readClasses(const char* classfile, double modthresh) {
  FILE* fp = fopen(classfile, "rb");
  if (fp == NULL) {
	  printf("Error while reading edge class file: could not open file %s\n", classfile); 
	  return false;
  }
  if (fread(&nclass, sizeof(int), 1, fp) != 1) {
	  fclose(fp);
	  printf("Error while reading edge class file: could not read 'nclass': int\n"); 
	  return false;
  }
  classes = (int*) malloc(2*m*sizeof(int));
  if (fread(classes, sizeof(int), 2*m, fp) != 2 * m) {
	  fclose(fp);
	  free(classes);
	  printf("Error while reading edge class file: could not read 'classes': int[2*edgeWeightsNumber]"); 
	  return false;
  }
  fclose(fp);
  
  if (modthresh >= 0) {

    // threshold for weights from when we treat it as modular
    double wthresh = 0.05 + 0.95*exp(-modthresh); 

    int ct=0;
    for (unsigned int i=0; i<m; i++) {
      if (weights[i] > wthresh) {
	classes[2*i] = nclass;
	classes[2*i+1] = nclass;
	++ct;
      }
    }
    printf("%d of %d not submod: %1.2e\n", ct, m, ((double)ct/(double)m)*100);

  }
  return true;
}


// set classes
// if you do not want modular treatment, then set modthresh < 0
// classfile should be of the form:
// num_classes
// int array of edge classes, length 2*m
bool Edat_read::readClasses(int* classes_in, int num_classes, double modthresh) {

  classes = classes_in;
  nclass = num_classes;
  
  if (modthresh >= 0) {

    // threshold for weights from when we treat it as modular
    double wthresh = 0.05 + 0.95*exp(-modthresh); 

    int ct=0;
    for (unsigned int i=0; i<m; i++) {
      if (weights[i] > wthresh) {
	classes[2*i] = nclass;
	classes[2*i+1] = nclass;
	++ct;
      }
    }
    printf("%d of %d not submod: %1.2e\n", ct, m, ((double)ct/(double)m)*100);

  }
  return false;
}
