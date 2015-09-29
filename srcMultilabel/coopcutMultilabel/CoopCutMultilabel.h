#pragma once

#include <vector>

#include "../gcoLibrary/gco-v3.0/GCoptimization.h"

class CoopCutMultilabel{
	public:
		CoopCutMultilabel();
		~CoopCutMultilabel();

		bool readUnaryTermsFromFile(const char* fileName);
		bool readPairwiseTermsFromFile(const char* fileName);

		void setLambda(double lambda) { lambda_ = lambda; }
		void setAlpha(double alpha) { alpha_ = alpha; }
		void setTheta(double theta);
		void setGreedyMaxIter(int maxIter) { greedyMaxIter_ = maxIter; }

		double getEnergy() const { return energy_; } 

		bool writeSolutionToFile(const char * fileName) const;

		bool mimimizeEnergy_greedy();
		bool minimizeEnergy_alphaExpansion();

	private:

		// model parameters
		double lambda_;
		double edgeClassTheta_;
		double alpha_;
		int greedyMaxIter_;

		std::vector<double> classThreshold_;

		// graph data
		int nodeNumber_;
		int edgeNumber_;
		int edgeClassNumber_;
		int labelNumber_;

		double* unaryTerms_;
		int getUnaryTermIndex(const int iNode, const int iLabel) const { return iLabel + iNode * labelNumber_; }
		
		int* edges_;
		int getEdgeEndIndex(const int iEdge, const int iDirection) const { return iEdge * 2 + iDirection; } 

		double* edgeWeights_;

		int* edgeClasses_;
		int getEdgeClassIndex(const int iEdge, const int iDirection) const {  return iEdge * 2 + iDirection; }
		
		// result data
		std::vector<int> labeling_;
		double energy_;

		// loading data
		bool unaryLoaded_;
		bool pairwiseLoaded_;
		bool resultsComputed_;
		double elapsedTime_;

		// internal functions
		int getFileSize(const char* filename) const;

		// callback for alpha-exp that ignore cooperation
		GCoptimization::EnergyTermType computePairwiseCost(
					GCoptimization::SiteID s1, GCoptimization::SiteID s2,
					GCoptimization::LabelID l1, GCoptimization::LabelID l2
					);
		friend	GCoptimization::EnergyTermType globalPairwiseCost(
					GCoptimization::SiteID s1, 	GCoptimization::SiteID s2, 
					GCoptimization::LabelID l1,	GCoptimization::LabelID l2,
					void * objectPtr
				);

		// callbacks that take coopreration into account
		GCoptimization::EnergyTermType computePairwiseCostOneBreakPoint(
					GCoptimization::SiteID s1, GCoptimization::SiteID s2,
					GCoptimization::LabelID l1, GCoptimization::LabelID l2
					);
		friend	GCoptimization::EnergyTermType globalPairwiseCostOneBreakPoint(
					GCoptimization::SiteID s1, 	GCoptimization::SiteID s2, 
					GCoptimization::LabelID l1,	GCoptimization::LabelID l2,
					void * objectPtr
				);
		
		void constructIndexingDataStructure(int maxNeighborNumber);
		int getEdgeIndex(int node0, int node1) const;
		
		int *neighborNumber_;
		int *allNeighboringNodes_;
		double *allNeighboringWeights_;
		int *allNeighboringEdgeGlobalIndex_;

		int **neighboringNode_;
		double **neighboringEdgeWeight_;
		int **neighboringEdgeGlobalIndex_;

		// functions for mimimizeEnergy_greedy
		std::vector<std::vector< int >> hValue_;
		double oneRunAlphaExpansion(std::vector<int> &labeling);

		double computeCoopEnergy(const std::vector<int> labeling);

};