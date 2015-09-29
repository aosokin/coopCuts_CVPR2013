//This class holds the data that defines coopcut and reads this data from files

#ifndef __COOPCUTDATA_H__
#define __COOPCUTDATA_H__

#include <vector>
using std::vector;

// QPBO header
#include "../QPBO/QPBO-v1.3.src/QPBO.h"

// GraphCut header
#include "../graphCuts/maxflow-v3.02.src/graph.h"
typedef Graph<double, double, double> GraphCutObject;


class CoopCutMinimizer {
	// this class is implemented only for binary energy
public:
	CoopCutMinimizer(int nodeNumber = 0, int edgeNumber = 0, int edgeClassNumber= 0);
	~CoopCutMinimizer();

	// IO functions for Stefanie's format
	bool readUnaryPotentials_Stefanie(const char* unaryFile);
	bool readPairwisePotentials_Stefanie(const char* edgeWeightFile, const char* edgeClassFile, bool readDiagonal = false);
	bool writeLabeling_Stefanie(const char* labelingFile, const double elapsedTime = 0.0) const;

	// Setting method parameters
	void setClassTheta(double edgeClassTheta, double edgeClassLastTheta = 1.0) { edgeClassTheta_ = edgeClassTheta; edgeClassLastTheta_ = edgeClassLastTheta; }
	void setLambda(double lambda) { lambda_ = lambda; }
	void setAlpha(double alpha) { alpha_ = alpha; }
	void setVerbosity(int level) { verbosityLevel_ = level;} // 0 - no printing; 1 - final printing; 2 - iter printing

	// Minimizing algorithms
	void minimize_QPBO();
	void minimize_modular_GraphCut();
	void minimize_exhaustive_GraphCut();
	void minimize_exhaustive_GraphCut_dynamic();
	void minimize_greedy_marginals(int maxIter = 10);
	void minimize_allsteps_marginals(int maxIter = 10);

	// Result getters
	double getEnergy() const { return energy_  + energyConstant_; }
	double getLowerBound() const { return lowerBound_  + energyConstant_; }
	int getNumberNodesUnlabeled() const { return numNodesUnlabeled_; }
	int getLabel(int iNode) const { return labeling_[iNode]; }


private:
	//parameters
	double edgeClassTheta_;
	double edgeClassLastTheta_;
	double lambda_; // weight for pairwise potentials
	double alpha_; // slope after the breakpoint
	int verbosityLevel_; // level of printing

	// qualitative charasteristics
	int nodeNumber_;
	double energyConstant_;
	int edgeClassNumber_;
	int edgeNumber_;

	// data
	double* unaryPotentials_;
	double* classTreshold_;
	double* totalEdgeClassWeight_;

	int* edges_;
	inline int getEdgeNodeIndex(int iEdge, int iNode) const { return 2 * iEdge + iNode; }
	
	int* edgeClass_;
	inline int getEdgeClassIndex(int iEdge, int direction) const { return 2 * iEdge + direction; }

	double* edgeWeights_;
	inline int getEdgeWeightIndex(int iEdge, int direction) const { return iEdge; }
	
	//Stefanie constant = 50.0
	const double constantStefanie_;

	// Method results
	bool resultsComputed_;
	double lowerBound_;
	double energy_;
	vector<int> labeling_;
	int numNodesUnlabeled_;

	// Computes the energy for the labeling
	// CAUTION: this version ignores the energyConstant_
	double computeEnergy(const vector<int>& labeling) const;	


private:
	int getFileSize(const char* filename) const;
	bool constructEdges_Stefanie(int height, int width);

	// functions to set up an order for exhaustive search in minimize_exhastive_GraphCut()
	void setOrderOfClassSwitches(const vector<double>& classPriority, vector<int>* switches) const; 
	void createHanoiTowerArray(int N, vector<int> *answer) const;

	// function to dynamically change the graph for dynamic maxflow
	double updateTheGraph(GraphCutObject* graphCut, int switchDirection, const vector<GraphCutObject::arc_id> &arcs, const vector<double> &weights) const;

	

};

typedef std::pair<double,int> PairDoubleInt;

inline bool less_PairDoubleInt ( const PairDoubleInt& l, const PairDoubleInt& r)
{
	return l.first < r.first; 
}

#endif