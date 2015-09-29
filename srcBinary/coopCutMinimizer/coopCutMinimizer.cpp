
#include "coopCutMinimizer.h"
#include <stdio.h>

#include <sys/stat.h>


CoopCutMinimizer::CoopCutMinimizer(int nodeNumber, int edgeNumber, int edgeClassNumber)
	: nodeNumber_(nodeNumber), 
	edgeNumber_(edgeNumber), 
	edgeClassNumber_(edgeClassNumber),
	unaryPotentials_(NULL),
	energyConstant_(0),
	classTreshold_(NULL),
	edges_(NULL),
	edgeClass_(NULL),
	edgeWeights_(NULL),
	resultsComputed_(false),
	lowerBound_(0),
	energy_(0),
	labeling_(nodeNumber_, -1),
	numNodesUnlabeled_(-1),
	totalEdgeClassWeight_(NULL),
	constantStefanie_(50.0)	,
	lambda_(1.0),
	edgeClassTheta_(0.5),
	edgeClassLastTheta_(1.0),
	alpha_(0),
	verbosityLevel_(0)
{
}

CoopCutMinimizer::~CoopCutMinimizer() {
	if (unaryPotentials_ != NULL) free(unaryPotentials_);
	if (classTreshold_ != NULL) free(classTreshold_);
	if (edges_ != NULL) free(edges_);
	if (edgeClass_ != NULL) free(edgeClass_);
	if (edgeWeights_ != NULL) free(edgeWeights_);
	if (totalEdgeClassWeight_ != NULL) free(totalEdgeClassWeight_);	
}

bool CoopCutMinimizer::readUnaryPotentials_Stefanie(const char* unaryFile) {
	//nodeNumber_ should be set beforehand
	if (nodeNumber_ <= 0) { return false;}

	// get the size of the file
	int fileSize = getFileSize(unaryFile);
	if (fileSize < 2 * sizeof(double) * nodeNumber_) {
		return false;
	}

	//open the file
	FILE* fp = fopen(unaryFile, "rb");
	if(fp == NULL) {return false;}

	//create the buffer
	double* unaryWeights = (double*) malloc(2 * nodeNumber_ * sizeof(double));
	
	//read the file ti the buffer
	size_t elementsRead = fread(unaryWeights, sizeof(double), 2 * nodeNumber_, fp); // all s-weights and then all t-weights
	fclose(fp);

	// check if everything was succefully read
	if (elementsRead != 2 * nodeNumber_) {	free(unaryWeights);	return false;	}
	
	// add unary potentials
	if (unaryPotentials_ != NULL) free(unaryPotentials_);
	unaryPotentials_ = (double*) malloc(nodeNumber_ * sizeof(double));

	for(int iNode = 0; iNode < nodeNumber_; ++iNode) {
		unaryPotentials_[iNode] = (unaryWeights[iNode + nodeNumber_] - unaryWeights[iNode]) / lambda_;
		energyConstant_ += unaryWeights[iNode] / lambda_;
	}
	free(unaryWeights);	
	unaryWeights = NULL;

	return true;
}


bool CoopCutMinimizer::readPairwisePotentials_Stefanie(const char* edgeWeightFile, const char* edgeClassFile, bool readDiagonal) {
	// get the size of the edgeWeightFile
	int edgeWeightFileSize = getFileSize(edgeWeightFile);
	if (edgeWeightFileSize < 3 * sizeof(int)) { return false; } // Check if file is extremely small
	
	// read the edge weights (created by the matlab code, makeweights.m)
	FILE* fp = fopen(edgeWeightFile,"rb");
	if(fp == NULL) {return false;}
	int height, width, weightsNumber;
	if (fread(&height, sizeof(int), 1, fp) != 1) { fclose(fp); return false; }
	if (fread(&width, sizeof(int), 1, fp) != 1) { fclose(fp); return false; }
	if (fread(&weightsNumber, sizeof(int), 1, fp) != 1) { fclose(fp); return false; }
	double* weights = (double*) malloc(weightsNumber * sizeof(double));
	if (fread(weights, sizeof(double), weightsNumber, fp) != weightsNumber) { fclose(fp); free(weights); return false; }
	fclose(fp);
	
	// set the number of nodes
	nodeNumber_ = height * width;

	// set the number of edges and t
	edgeNumber_ = (height - 1) * width + height * (width - 1);
	if (readDiagonal) {
		edgeNumber_ += 2 * (height - 1) * (width - 1);
	}

	// check if file is large enough
	if (weightsNumber < edgeNumber_) { free(weights); return false; }

	// add the edge weights: symmetric
	if (edgeWeights_ != NULL) free(edgeWeights_);
	edgeWeights_ = weights;
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		edgeWeights_[iEdge] = weights[iEdge] * constantStefanie_ ;
	}
	weights = NULL;

	// get the size of the edgeClassFile
	int edgeClassFileSize = getFileSize(edgeClassFile);
	if ( edgeClassFileSize < sizeof(int) + 2 * sizeof(int) * edgeNumber_) {	return false;} //Check if file is long enough

	// read the edge classes
	fp = fopen(edgeClassFile,"rb");
	if(fp == NULL) {return false;}
	if (fread(&edgeClassNumber_, sizeof(int), 1, fp) != 1) { fclose(fp); return false; }
	int* edgeClasses = (int*) malloc(edgeNumber_ * 2 * sizeof(int));

	size_t numberOfReadEdgeClasses = fread(edgeClasses, sizeof(int), edgeNumber_ * 2, fp);
	if (numberOfReadEdgeClasses != edgeNumber_ * 2) { fclose(fp); free(edgeClasses); return false; }
	fclose(fp);

	// add modular class to the  class set
	++edgeClassNumber_;

	// analize the edge classes
	
	// create arrays for edge thresholds
	if (classTreshold_ != NULL) free(classTreshold_);
	classTreshold_ = (double*) malloc(edgeClassNumber_ * sizeof(double));
	if (totalEdgeClassWeight_ != NULL) free(totalEdgeClassWeight_);
	totalEdgeClassWeight_ = (double*) malloc(edgeClassNumber_ * sizeof(double));
	//initialize totalEdgeClassWeight_
	for(int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		totalEdgeClassWeight_[iClass] = 0;
	
	// sum up the weights of the classes
	// save classes of edges
	if( edgeClass_ != NULL ) free(edgeClass_);
	edgeClass_ = (int*) malloc(2 * edgeNumber_ * sizeof(int));

	for( int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		double curThresh;
	
		//add forward edge and save its class
		int curClass = (edgeClass_[getEdgeClassIndex(iEdge, 0)] = edgeClasses[2 * iEdge + 1]);

		// adjust thereshhold for a forward edge
		totalEdgeClassWeight_[ curClass ] += edgeWeights_[getEdgeWeightIndex(iEdge, 0)];
		
		//add backward edge
		curClass = (edgeClass_[getEdgeClassIndex(iEdge, 1)] = edgeClasses[2 * iEdge]);
		
		// adjust thereshhold for a backward edge
		totalEdgeClassWeight_[ curClass ] += edgeWeights_[getEdgeWeightIndex(iEdge, 1)];
	}
	free(edgeClasses);
	edgeClasses = NULL;

	// set up the class thresholds
	for(int iClass = 0; iClass + 1 < edgeClassNumber_; ++iClass) {
		classTreshold_[iClass] = totalEdgeClassWeight_[iClass] * edgeClassTheta_;
	}
	classTreshold_[edgeClassNumber_ - 1] = totalEdgeClassWeight_[edgeClassNumber_ - 1] * edgeClassLastTheta_;

	// connect edges with nodes: depends on format of data very much
	return constructEdges_Stefanie(height, width);
}


int CoopCutMinimizer::getFileSize(const char* filename) const
{
    struct stat stat_buf;
    int rc = stat(filename, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

bool CoopCutMinimizer::constructEdges_Stefanie(int height, int width)
{
	// reserve memory for the edge array
	if(edges_ != NULL) free(edges_);
	edges_ = (int*) malloc(2 * edgeNumber_ * sizeof(int));

	// add the edges
	int edgeCount = 0;

	int curBlockEdgeNumber = width * (height - 1);

	//add vertical edges
	if (edgeNumber_ >= curBlockEdgeNumber) {
		for (int iCol = 0; iCol < width; ++iCol) 
			for (int iRow = 0; iRow < height - 1; ++iRow) {
				edges_[getEdgeNodeIndex(edgeCount, 0)] = iRow * width + iCol;
				edges_[getEdgeNodeIndex(edgeCount, 1)] = (iRow + 1) * width + iCol;
				++edgeCount;
			}

		if (edgeNumber_ == edgeCount) {
			return true;
		}
	}

	//add horizontal edges
	curBlockEdgeNumber += (width - 1) * height;
	if (edgeNumber_ >= curBlockEdgeNumber) {
		for (int iCol = 0; iCol < width - 1; ++iCol) 
			for (int iRow = 0; iRow < height; ++iRow) {
				edges_[getEdgeNodeIndex(edgeCount, 0)] = iRow * width + iCol + 1;
				edges_[getEdgeNodeIndex(edgeCount, 1)] = iRow * width + iCol;
				++edgeCount;
			}

		if (edgeNumber_ == edgeCount) {
			return true;
		}
	}

	//add diagonal '\' edges
	curBlockEdgeNumber += (width - 1) * (height - 1);
	if (edgeNumber_ >= curBlockEdgeNumber) {
		for (int iCol = 0; iCol < width - 1; ++iCol) 
			for (int iRow = 0; iRow < height - 1; ++iRow) {
				edges_[getEdgeNodeIndex(edgeCount, 0)] = iRow * width + iCol ;
				edges_[getEdgeNodeIndex(edgeCount, 1)] = (iRow + 1) * width + iCol + 1;
				++edgeCount;
			}

		if (edgeNumber_ == edgeCount) {
			return true;
		}
	}

	//add diagonal '/' edges
	curBlockEdgeNumber += (width - 1) * (height - 1);
	if (edgeNumber_ >= curBlockEdgeNumber) {
		for (int iCol = 0; iCol < width - 1; ++iCol) 
			for (int iRow = 1; iRow < height; ++iRow) {
				edges_[getEdgeNodeIndex(edgeCount, 0)] = iRow * width + iCol ;
				edges_[getEdgeNodeIndex(edgeCount, 1)] = (iRow - 1) * width + iCol + 1;
				++edgeCount;
			}
		if (edgeNumber_ == edgeCount) {
			return true;
		}
	}
	return false;
}

bool CoopCutMinimizer::writeLabeling_Stefanie(const char* labelingFile, const double elapsedTime) const
{
	//check if there is a solution
	if (!resultsComputed_) return false;

	//create the buffer array
	signed char* buffer = (signed char*) malloc(nodeNumber_ * sizeof(signed char));
	for(int iNode = 0; iNode < nodeNumber_; ++iNode)
		buffer[iNode] = (signed char)labeling_[iNode];

	// write out node labels
	FILE* of = fopen(labelingFile, "wb");
	fwrite(&nodeNumber_, sizeof(int), 1, of);
	fwrite(buffer, sizeof(signed char), nodeNumber_, of);

	// write elapsedTime to the file
	fwrite(&elapsedTime, sizeof(double), 1, of);

	// write energy to a file
	double energy = getEnergy();
	fwrite(&energy, sizeof(double), 1, of);

	fclose(of);
	free(buffer);
  
	return true;
}