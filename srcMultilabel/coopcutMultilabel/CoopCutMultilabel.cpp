
#include "CoopCutMultilabel.h"
#include "../gcoLibrary/gco-v3.0/GCoptimization.h"

#include <vector>
using std::vector;

#include <stdio.h>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <assert.h>

CoopCutMultilabel::CoopCutMultilabel() :
		nodeNumber_(-1),
		edgeNumber_(-1),
		edgeClassNumber_(-1),
		labelNumber_(-1),
		unaryTerms_(NULL),
		edges_(NULL),
		edgeWeights_(NULL),
		edgeClasses_(NULL),
		labeling_(0),
		unaryLoaded_(false),
		pairwiseLoaded_(false),
		resultsComputed_(false),
		elapsedTime_(-1),
		neighborNumber_(NULL),
		neighboringNode_(NULL),
		neighboringEdgeWeight_(NULL),
		allNeighboringNodes_(NULL),
		allNeighboringWeights_(NULL),
		neighboringEdgeGlobalIndex_(NULL),
		allNeighboringEdgeGlobalIndex_(NULL),
		lambda_( 1.0 ),
		edgeClassTheta_ ( 0.5 ),
		alpha_(0.1),
		classThreshold_(0),
		hValue_(0),
		greedyMaxIter_(10)
{
}

CoopCutMultilabel::~CoopCutMultilabel()
{
	if ( unaryTerms_ != NULL) delete [] unaryTerms_;
	if ( edges_ != NULL) delete [] edges_;
	if ( edgeWeights_ != NULL) delete [] edgeWeights_;
	if ( edgeClasses_ != NULL) delete [] edgeClasses_;

	if ( neighborNumber_ != NULL) delete [] neighborNumber_;
	if ( allNeighboringNodes_ != NULL) delete [] allNeighboringNodes_;
	if ( allNeighboringWeights_ != NULL) delete [] allNeighboringWeights_;
	if ( allNeighboringEdgeGlobalIndex_ != NULL) delete [] allNeighboringEdgeGlobalIndex_;
	if ( neighboringNode_ != NULL) delete [] neighboringNode_;
	if ( neighboringEdgeWeight_ != NULL) delete [] neighboringEdgeWeight_;
	if ( neighboringEdgeGlobalIndex_ != NULL) delete [] neighboringEdgeGlobalIndex_;

}

void CoopCutMultilabel::setTheta(double theta)
{
	edgeClassTheta_ = theta;

	if (pairwiseLoaded_) {
		for (int iClass = 0; iClass < edgeClassNumber_ + 1; ++iClass)
			classThreshold_[ iClass ] = 0;

		for (int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
			classThreshold_[ getEdgeClassIndex( iEdge, 0 ) ] += edgeWeights_[ iEdge ];
			classThreshold_[ getEdgeClassIndex( iEdge, 1 ) ] += edgeWeights_[ iEdge ];
		}

		for (int iClass = 0; iClass < edgeClassNumber_; ++iClass)
			classThreshold_[ iClass ] *= edgeClassTheta_;
	}
}

int CoopCutMultilabel::getFileSize(const char* filename) const
{
    struct stat stat_buf;
    int rc = stat(filename, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}


bool CoopCutMultilabel::readUnaryTermsFromFile(const char* fileName)
{
	// get the size of the file
	int fileSize = getFileSize(fileName);
	if (fileSize < 2 * sizeof(int)) {
		return false;
	}

	//open the file
	FILE* fp = fopen(fileName, "rb");
	if(fp == NULL) {return false;}

	// read the header
	fread(&nodeNumber_, sizeof(int), 1, fp);
	fread(&labelNumber_, sizeof(int), 1, fp);

	if (labelNumber_ <= 0 || nodeNumber_ <= 0) {
		fclose( fp );
		return false;
	}

	//create the buffer
	double* unaryWeights = (double*) malloc(labelNumber_ * nodeNumber_ * sizeof(double));
	
	//read the file ti the buffer
	size_t elementsRead = fread(unaryWeights, sizeof(double), labelNumber_ * nodeNumber_, fp); // all unaryWeights
	fclose(fp);

	// check if everything was succefully read
	if (elementsRead != labelNumber_ * nodeNumber_) {	free(unaryWeights);	return false;	}
	
	// add unary potentials
	unaryTerms_ = unaryWeights;
	unaryWeights = NULL;

	unaryLoaded_ = true;
	return true;
}

bool CoopCutMultilabel::readPairwiseTermsFromFile(const char* fileName)
{
	// get the size of the file
	int fileSize = getFileSize(fileName);
	if (fileSize < 2 * sizeof(int)) {
		return false;
	}

	//open the file
	FILE* fp = fopen(fileName, "rb");
	if(fp == NULL) {return false;}

	// read the header
	fread(&edgeNumber_, sizeof(int), 1, fp);
	fread(&edgeClassNumber_, sizeof(int), 1, fp);

	if (edgeNumber_ <= 0 || edgeClassNumber_ <= 0) {
		fclose( fp );
		return false;
	}

	// initialize counters for the total sums of weights of the current class
	classThreshold_.clear();
	classThreshold_.resize(edgeClassNumber_ + 1, 0);

	// go though all edges
	int iNode[2];
	int edgeClass[2];
	double edgeWeight;
	size_t elementsRead;

	edges_ = (int*) malloc( 2 * edgeNumber_ * sizeof(int) );
	edgeWeights_ = (double*) malloc( edgeNumber_ * sizeof(double) );
	edgeClasses_  = (int*) malloc( 2 * edgeNumber_ * sizeof(int) );
		
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		//read the file to the buffer
		elementsRead = fread(iNode, sizeof(int), 2, fp); // all edge ends
		if (elementsRead != 2) { fclose(fp); return false; }
		
		elementsRead = fread(&edgeWeight, sizeof(double), 1, fp); // edge weight
		if (elementsRead != 1) { fclose(fp); return false; }
		
		elementsRead = fread(edgeClass, sizeof(int), 2, fp); // edge class
		if (elementsRead != 2) { fclose(fp); return false; }
		
		edges_[ getEdgeEndIndex(iEdge, 0) ] =  iNode[ 0 ];
		edges_[ getEdgeEndIndex(iEdge, 1) ] =  iNode[ 1 ];

		edgeWeights_[ iEdge ] = edgeWeight;
		
		edgeClasses_[ getEdgeClassIndex(iEdge, 0) ] = edgeClass[ 0 ];
		edgeClasses_[ getEdgeClassIndex(iEdge, 1) ] = edgeClass[ 1 ];

		classThreshold_[ edgeClass[ 0 ] ] += edgeWeight;
		classThreshold_[ edgeClass[ 1 ] ] += edgeWeight;
	}

	for(int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		classThreshold_[iClass] *= edgeClassTheta_;	

	pairwiseLoaded_ = true;
	return true;
}

		
bool CoopCutMultilabel::writeSolutionToFile(const char * fileName) const
{
	//check if there is a solution
	if (!resultsComputed_) return false;

	//create the buffer array
	signed char* buffer = (signed char*) malloc(nodeNumber_ * sizeof(signed char));
	for(int iNode = 0; iNode < nodeNumber_; ++iNode)
		buffer[iNode] = (signed char)labeling_[iNode];

	// write out node labels
	FILE* of = fopen(fileName, "wb");
	fwrite(&nodeNumber_, sizeof(int), 1, of);
	fwrite(buffer, sizeof(signed char), nodeNumber_, of);

	// write elapsedTime to the file
	fwrite(&elapsedTime_, sizeof(double), 1, of);

	// write energy to a file
	double energy = getEnergy();
	fwrite(&energy, sizeof(double), 1, of);

	fclose(of);
	free(buffer);
  
	return true;
}

GCoptimization::EnergyTermType globalPairwiseCost
	(
	GCoptimization::SiteID s1, 
	GCoptimization::SiteID s2, 
	GCoptimization::LabelID l1, 
	GCoptimization::LabelID l2,
	void * objectPtr
	)
{
	CoopCutMultilabel* objectPtrCorType = (CoopCutMultilabel*) objectPtr;
	return objectPtrCorType -> computePairwiseCost(s1, s2, l1, l2);
}

GCoptimization::EnergyTermType CoopCutMultilabel::computePairwiseCost
	(
	GCoptimization::SiteID s1, 
	GCoptimization::SiteID s2, 
	GCoptimization::LabelID l1, 
	GCoptimization::LabelID l2
	)
{
	if (l1 == l2) return 0;
	
	int edgeForward = getEdgeIndex(s1, s2); 

	assert(neighboringEdgeGlobalIndex_ != NULL);

	int edgeGlobalIndex = neighboringEdgeGlobalIndex_[s1][edgeForward];
	
	return edgeWeights_[ edgeGlobalIndex ] * lambda_;
}

GCoptimization::EnergyTermType globalPairwiseCostOneBreakPoint
	(
	GCoptimization::SiteID s1, 
	GCoptimization::SiteID s2, 
	GCoptimization::LabelID l1, 
	GCoptimization::LabelID l2,
	void * objectPtr
	)
{
	CoopCutMultilabel* objectPtrCorType = (CoopCutMultilabel*) objectPtr;
	return objectPtrCorType -> computePairwiseCostOneBreakPoint(s1, s2, l1, l2);
}

GCoptimization::EnergyTermType CoopCutMultilabel::computePairwiseCostOneBreakPoint
	(
	GCoptimization::SiteID s1, 
	GCoptimization::SiteID s2, 
	GCoptimization::LabelID l1, 
	GCoptimization::LabelID l2
	)
{
	if (l1 == l2) return 0;
	
	int edgeForward = getEdgeIndex(s1, s2); 

	assert(neighboringEdgeGlobalIndex_ != NULL);
	int edgeGlobalIndex = neighboringEdgeGlobalIndex_[s1][edgeForward];

	assert( edgeGlobalIndex  == neighboringEdgeGlobalIndex_[s2][getEdgeIndex(s2, s1)] );

	if (s1 != edges_[ getEdgeEndIndex(edgeGlobalIndex, 0) ] ){
		int c = s1;
		s1 = s2;
		s2 = c;

		c = l1;
		l1 = l2;
		l2 = c;
	}

	assert(s1 == edges_[ getEdgeEndIndex(edgeGlobalIndex, 0) ]);
	assert(s2 == edges_[ getEdgeEndIndex(edgeGlobalIndex, 1) ]);



	double edgeWeight = edgeWeights_[ edgeGlobalIndex ];
	int edgeClassForward = edgeClasses_[ getEdgeClassIndex(edgeGlobalIndex, 0) ];
	int edgeClassBackward = edgeClasses_[ getEdgeClassIndex(edgeGlobalIndex, 1) ];

	double weightForward = 0;
	double weightBackward = 0;

	if (edgeClassForward == edgeClassNumber_) {
		weightForward = edgeWeight;
	} else {
		weightForward = (hValue_[ l1 ][ edgeClassForward ] == 1) ? edgeWeight : edgeWeight * alpha_;
	}

	if (edgeClassBackward == edgeClassNumber_) {
		weightBackward = edgeWeight;
	} else {
		weightBackward = (hValue_[ l2 ][ edgeClassBackward ] == 1) ? edgeWeight : edgeWeight * alpha_;
	}

	return lambda_ * (weightForward + weightBackward) / 2;
}


void CoopCutMultilabel::constructIndexingDataStructure( int maxNeighborNumber )
{
	neighborNumber_ = new int[ nodeNumber_ ];
	allNeighboringNodes_ = new int[nodeNumber_ * maxNeighborNumber];
	allNeighboringWeights_ = new double[nodeNumber_ * maxNeighborNumber];
	allNeighboringEdgeGlobalIndex_ = new int[nodeNumber_ * maxNeighborNumber];
	
	neighboringNode_ = new int*[nodeNumber_];
	neighboringEdgeWeight_ = new double*[nodeNumber_];
	neighboringEdgeGlobalIndex_ = new int*[nodeNumber_];

	for(int iNode = 0; iNode < nodeNumber_; ++iNode){
		neighborNumber_[ iNode ] = 0;
		neighboringNode_[ iNode ] = allNeighboringNodes_ + (iNode * maxNeighborNumber);
		neighboringEdgeWeight_[ iNode ] = allNeighboringWeights_ + (iNode * maxNeighborNumber);
		neighboringEdgeGlobalIndex_[ iNode ] = allNeighboringEdgeGlobalIndex_ + (iNode * maxNeighborNumber);
	}

	// go through all edges and add them
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		//	get all data
		int node0 = edges_[ getEdgeEndIndex(iEdge, 0) ];
		assert(node0 >= 0 && node0 < nodeNumber_);

		int node1 = edges_[ getEdgeEndIndex(iEdge, 1) ];
		assert(node1 >= 0 && node1 < nodeNumber_);

		double curWeight = 1.0;
		
		assert(neighborNumber_[node0] < maxNeighborNumber);
		assert(neighborNumber_[node1] < maxNeighborNumber);

		neighboringNode_[node0][ neighborNumber_[node0] ] = node1;
		neighboringNode_[node1][ neighborNumber_[node1] ] = node0;

		neighboringEdgeWeight_[node0][ neighborNumber_[node0] ] = curWeight;
		neighboringEdgeWeight_[node1][ neighborNumber_[node1] ] = curWeight;

		neighboringEdgeGlobalIndex_[node0][ neighborNumber_[node0] ] = iEdge;
		neighboringEdgeGlobalIndex_[node1][ neighborNumber_[node1] ] = iEdge;

		++neighborNumber_[node0];
		++neighborNumber_[node1];
	}
}

int CoopCutMultilabel::getEdgeIndex(int node0, int node1) const
{
	assert( neighborNumber_ != NULL );
	assert( neighboringNode_ != NULL );
	int i = 0;
	while (i < neighborNumber_[ node0 ] && neighboringNode_[ node0 ][i] != node1 ) 
		++i;

	assert( i < neighborNumber_[ node0 ]);

	return i;
}		


bool CoopCutMultilabel::minimizeEnergy_alphaExpansion()
{
	double tAlgoStart = clock();
	double tStart = clock();
	printf("Preparing energy\n");

	if ( !pairwiseLoaded_ || !unaryLoaded_ ) {
		return false;
	}

	labeling_.resize( nodeNumber_ );   // stores result of optimization

	// unaryTerms_  - array of unary potentials

	// prepare grid data structure
	const int maxNeighborNumber = 8; // this number is a hack!!! TODO: compute it from the input
	constructIndexingDataStructure( maxNeighborNumber ); 

	try{
	
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(nodeNumber_, labelNumber_);
		gc->setDataCost(unaryTerms_);
		
		gc->setSmoothCost( globalPairwiseCost, this );

		//gc->setSmoothCost(smooth);
		
		printf("Time: %f\n", (clock() - tStart) / CLOCKS_PER_SEC);

		// now set up a grid neighborhood system
		printf("Setting the neighborhood\n");
		tStart = clock();
		gc -> setAllNeighbors( neighborNumber_, neighboringNode_, neighboringEdgeWeight_);

		printf("Time: %f\n", (clock() - tStart) / CLOCKS_PER_SEC);
		
		printf("Before optimization energy is %f\n", gc->compute_energy());
		
		printf("Running alpha-expansion\n");
		tStart = clock();
		gc->expansion(10);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		printf("Time: %f\n", (clock() - tStart) / CLOCKS_PER_SEC);

		printf("Computing the energy and labeling\n");
		tStart = clock();

		energy_ = gc->compute_energy();
		printf("After optimization energy is %f\n", energy_);

		for ( int  iNode = 0; iNode < nodeNumber_; ++iNode )
			labeling_[ iNode ] = gc->whatLabel( iNode );
		printf("Time: %f\n", (clock() - tStart) / CLOCKS_PER_SEC);

		printf("Deleting GCO object\n");
		tStart = clock();

		delete gc;

		printf("Time: %f\n", (clock() - tStart) / CLOCKS_PER_SEC);

	}
	catch (GCException e){
		e.Report();
	}

	elapsedTime_ = (clock() - tAlgoStart) / CLOCKS_PER_SEC;
	resultsComputed_ = true;

	double energy = computeCoopEnergy( labeling_ );
	if( true || fabs( energy - energy_ ) > 1e-3 ){
		printf("Correct energy: %f\n", energy);
		energy_ = energy;
	}

	return true;
}

bool CoopCutMultilabel::mimimizeEnergy_greedy()
{
	double tAlgoStart = clock();
	double tStart = clock();
	printf("Preparing energy\n");

	if ( !pairwiseLoaded_ || !unaryLoaded_ ) {
		return false;
	}

	labeling_.resize( nodeNumber_ );   // stores result of optimization

	// unaryTerms_  - array of unary potentials

	// prepare grid data structure
	const int maxNeighborNumber = 8; // this number is a hack!!! TODO: compute it from the input
	constructIndexingDataStructure( maxNeighborNumber ); 
	

	// initialize hValue_
	hValue_.resize( labelNumber_, vector<int>( edgeClassNumber_, 1) );
	double hSum = 0.0;

	energy_ = oneRunAlphaExpansion( labeling_ );
	printf("Initial energy is %f\n", energy_);

	for(int iIter = 0; iIter < greedyMaxIter_; ++iIter) {
		//go through all h
		bool changed = false;

		for(int iLabel = 0; iLabel < labelNumber_; ++iLabel)
			for(int iEdgeClass = 0; iEdgeClass < edgeClassNumber_; ++ iEdgeClass)
			{
				// try to switch the current h
				hValue_[iLabel][iEdgeClass] = 1 - hValue_[iLabel][iEdgeClass];
				if (hValue_[iLabel][iEdgeClass] == 1) {
					hSum -= (1 - alpha_) * classThreshold_[iEdgeClass];
				} else {
					hSum += (1 - alpha_) * classThreshold_[iEdgeClass];
				}

				vector<int> curLabeling(0);

				double energy = oneRunAlphaExpansion(curLabeling);

				if (energy + lambda_ * hSum < energy_) {
					// accept the change 
					energy_ = energy + lambda_ * hSum;
					labeling_ = curLabeling;
					changed = true;

					printf("Changing h(%d, %d) to %d, energy is %f\n", iLabel, iEdgeClass, hValue_[iLabel][iEdgeClass], energy_);
				} else {
					// return h back

					printf("Changing h(%d, %d) to %d REJECTED, energy is %f\n", iLabel, iEdgeClass, hValue_[iLabel][iEdgeClass], energy + lambda_  * hSum);

					hValue_[iLabel][iEdgeClass] = 1 - hValue_[iLabel][iEdgeClass];
					if (hValue_[iLabel][iEdgeClass] == 1) {
						hSum -= (1 - alpha_) * classThreshold_[iEdgeClass];
					} else {
						hSum += (1 - alpha_) * classThreshold_[iEdgeClass];
					}
				}

			}

		if (!changed)
			break;
	}
	
	elapsedTime_ = (clock() - tAlgoStart) / CLOCKS_PER_SEC;
	printf("Energy: %f, time: %f\n", energy_, elapsedTime_);
	

	double energy = computeCoopEnergy( labeling_ );
	if( true || fabs( energy - energy_ ) > 1e-3 ){
		printf("Correct energy: %f\n", energy);
		energy_ = energy;
	}

	
	resultsComputed_ = true;

	return true;
}

double CoopCutMultilabel::oneRunAlphaExpansion(std::vector<int> &labeling) 
{
	double energy = 0.0;
	try{
	
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(nodeNumber_, labelNumber_);
		gc->setDataCost(unaryTerms_);
		gc->setSmoothCost( globalPairwiseCostOneBreakPoint, this );
		gc -> setAllNeighbors( neighborNumber_, neighboringNode_, neighboringEdgeWeight_);

		for(int iNode = 0; iNode < nodeNumber_; ++iNode)
			gc -> setLabel(iNode, 0);

		gc->expansion(10);
	
		energy = gc->compute_energy();
	
		labeling.resize(nodeNumber_);
		for ( int  iNode = 0; iNode < nodeNumber_; ++iNode )
			labeling[ iNode ] = gc->whatLabel( iNode );
		
		delete gc;
		
	}
	catch (GCException e){
		e.Report();
		return 1e+20;
	}

	return energy;
}

double CoopCutMultilabel::computeCoopEnergy(const std::vector<int> labeling)
{
	double energy = 0;

	for (int iNode = 0; iNode < nodeNumber_; ++iNode)
		energy += unaryTerms_[ getUnaryTermIndex(iNode, labeling[iNode]) ];

	vector<vector<double>> classCost(labelNumber_, vector<double>(edgeClassNumber_ + 1, 0) );
		
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		int class0 = edgeClasses_[ getEdgeClassIndex(iEdge, 0) ];
		int class1 = edgeClasses_[ getEdgeClassIndex(iEdge, 1) ];
		int node0 = edges_[ getEdgeEndIndex(iEdge, 0) ];
		int node1 = edges_[ getEdgeEndIndex(iEdge, 1) ];
		double curWeight = edgeWeights_[ iEdge ];

		if (labeling[node0] != labeling[node1]) {
			classCost[labeling[node0]][ class0 ] += curWeight / 2;
			classCost[labeling[node1]][ class1 ] += curWeight / 2;
		}
	}

	double pairwise = 0;
	for(int iLabel = 0; iLabel < labelNumber_; ++iLabel)
		for(int iEdgeClass = 0; iEdgeClass < edgeClassNumber_ + 1; ++ iEdgeClass)
		{
			double b = (1 - alpha_) * classThreshold_[iEdgeClass];

			double curW = classCost[iLabel][iEdgeClass];
			
			pairwise += std::min( curW, curW * alpha_ + b);
		}

	
	energy += lambda_ * pairwise;
	return energy;
}