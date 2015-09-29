#include "coopCutMinimizer.h"

#include <algorithm>
#include <stdio.h>
#include <time.h>

void CoopCutMinimizer::minimize_QPBO()
{
	QPBO<double>* minimizer = new QPBO<double>(nodeNumber_ + edgeClassNumber_ + 2 * edgeNumber_, 2 * 5 * edgeNumber_);

	// Add nodes: x_i and h_g
	minimizer -> AddNode(nodeNumber_ + edgeClassNumber_ + 2 * edgeNumber_); 
  
	// add unary potentials for x_i
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		minimizer -> AddUnaryTerm(iNode, 0,  unaryPotentials_[iNode]);
	}

	// add unary potentials for group nodes: h_g
	for (int iClass = 0; iClass < edgeClassNumber_; ++iClass){
	    minimizer -> AddUnaryTerm(nodeNumber_ + iClass, classTreshold_[iClass] * (1 - alpha_), 0); 
	}

	//Add edges
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		for(int iDirection = 0; iDirection <= 1; ++iDirection) {

			//find the weight of the current edge in the current direction
			double curWeight = edgeWeights_[getEdgeWeightIndex(iEdge, iDirection)];

			// find class of the current edge
			int curClass = edgeClass_[getEdgeClassIndex(iEdge, iDirection)];

			// find the indices of all nodes related to the current edge
			int x_i = edges_[getEdgeNodeIndex(iEdge, iDirection)];
			int x_j = edges_[getEdgeNodeIndex(iEdge, 1 - iDirection)];
			int h_g = nodeNumber_ + curClass;
			int z_ij = nodeNumber_ + edgeClassNumber_ + 2 * iEdge + iDirection;



			// if current class is modular do not add new edges
			if (classTreshold_[curClass] < (1.0 - (1e-2)) * totalEdgeClassWeight_[curClass]) {
				// unary for z_ij
				minimizer -> AddUnaryTerm(z_ij, 0, 2 * curWeight * ( 1.0 - alpha_));

				// edge z_ij to h_g
				minimizer -> AddPairwiseTerm(z_ij, h_g,	0, 0, 0, -curWeight * ( 1.0 - alpha_)); 
			
				// edge z_ij to x_i
				minimizer -> AddPairwiseTerm(z_ij, x_i, 0, 0, 0, -curWeight * ( 1.0 - alpha_));  
			
				// edge z_ij to x_j
				minimizer -> AddPairwiseTerm(z_ij, x_j,	0, 0, 0, -curWeight * ( 1.0 - alpha_));  
			
				// edges h_g to x_i -- positive
				minimizer -> AddPairwiseTerm(h_g, x_i, 0, 0, 0, curWeight * ( 1.0 - alpha_));

				// edge x_i to x_j
				minimizer -> AddPairwiseTerm(x_i, x_j, 0, 0, curWeight * alpha_, 0);
			} else {
				// add edge x_i to x_j
				minimizer -> AddPairwiseTerm(x_i, x_j, 0, 0, curWeight, 0);
			}
		}
	}

	// Merge parallel edges
	minimizer -> MergeParallelEdges();

	double tStart = clock();
	// Run the QPBO solver
	minimizer -> Solve();

	double elapsedTime = (clock() - tStart) / CLOCKS_PER_SEC;
	if (verbosityLevel_ >= 1) 
		printf("Time for QPBO run: %fs\n", elapsedTime);

	// Compute the weak persistent solution
	minimizer -> ComputeWeakPersistencies();

	// Get the lower bound
	lowerBound_ = minimizer -> ComputeTwiceLowerBound() / 2.0;

	// Get the energy filled with zeros
	energy_ = minimizer -> ComputeTwiceEnergy() / 2.0;
	
	// Get the solution and compute the number of unlabeled nodes
	numNodesUnlabeled_ = 0;
	labeling_.resize(nodeNumber_, -1);
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		labeling_[iNode] = minimizer -> GetLabel(iNode);
		if (labeling_[iNode] < 0) {
			++numNodesUnlabeled_;
		}
	}

	resultsComputed_ = true;
	
	delete minimizer;
}


void CoopCutMinimizer::minimize_modular_GraphCut()
{
	GraphCutObject* graphCut = new GraphCutObject(nodeNumber_, edgeNumber_);

	// Add nodes
	graphCut -> add_node(nodeNumber_);
  
	// add unary potentials for x_i
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		graphCut -> add_tweights(iNode, unaryPotentials_[iNode], 0);
	}

	//Add edges
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		//find the weight of the current edge in the current direction
		double curWeight0 = edgeWeights_[getEdgeWeightIndex(iEdge, 0)];
		double curWeight1 = edgeWeights_[getEdgeWeightIndex(iEdge, 1)];

		// find the indices of all nodes related to the current edge
		int x_i = edges_[getEdgeNodeIndex(iEdge, 0)];
		int x_j = edges_[getEdgeNodeIndex(iEdge, 1)];

		// Add pairwise terms
		graphCut -> add_edge(x_i, x_j, curWeight0, curWeight1);
	}
	//compute time only for maxFlow
	double flowTime = clock();
	// Run the GC solver
	double energyMaxflow = graphCut -> maxflow();
	
	flowTime = ( clock() - flowTime ) / CLOCKS_PER_SEC; 

	if (verbosityLevel_ >= 1) 
		printf("Flow time: %fs\n", flowTime);


	// Get the solution and compute the number of unlabeled nodes
	labeling_.resize(nodeNumber_, -1);
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		labeling_[iNode] = graphCut -> what_segment(iNode);
	}
	
	// The current energy is ignoring the cooperation potentials: thus recompute the energy
	energy_ = computeEnergy(labeling_);

	resultsComputed_ = true;
	numNodesUnlabeled_ = 0;
	
	delete graphCut;
}

void CoopCutMinimizer::minimize_exhaustive_GraphCut()
{
	//Determine priority of classes for an exhaustive search: -1.0 - not included, largest number - the most expensive change
	vector<double> classPriority(edgeClassNumber_, -1.0);
	for (int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		if (classTreshold_[iClass] < (1.0 - (1e-2)) * totalEdgeClassWeight_[iClass]) { // class is not modular
			classPriority[iClass] = totalEdgeClassWeight_[iClass];
		}

	// create a switch order for an exhaustive search
	vector<int> orderOfClassSwiches(0, 0);
	setOrderOfClassSwitches(classPriority, &orderOfClassSwiches);

	// inialize the starting switch configuration: each switch corresponds to h_g
	vector<int> classSwitches(edgeClassNumber_, 1); // Class that are modular are always on

	// compute the sum of current values of h_g
	double hSum = 0;
	for(int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		hSum += (1 - classSwitches[iClass]) * classTreshold_[iClass]  * (1 - alpha_);

	//create a graphCut object
	GraphCutObject* graphCut = new GraphCutObject(nodeNumber_, edgeNumber_);

	//set energy to infinity
	energy_ = 1e+40;

	// set the timer for flow computetions:
	double flowTime = 0, tic;

	
	for(int iIteration = 0; iIteration <= orderOfClassSwiches.size(); ++iIteration) {
		// Add nodes
		graphCut -> add_node(nodeNumber_);
  
		// add unary potentials for x_i
		for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
			graphCut -> add_tweights(iNode, unaryPotentials_[iNode], 0);
		}

		//Add edges
		for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
			//find the classes of current edges
			int curClass0 = edgeClass_[getEdgeClassIndex(iEdge, 0)];
			int curClass1 = edgeClass_[getEdgeClassIndex(iEdge, 1)];

			//find the weight of the current edge in the current direction
			double curWeight0 = (classSwitches[curClass0] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 0)] : alpha_ * edgeWeights_[getEdgeWeightIndex(iEdge, 0)];
			double curWeight1 = (classSwitches[curClass1] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 1)] : alpha_ * edgeWeights_[getEdgeWeightIndex(iEdge, 1)];
		
			// find the indices of all nodes related to the current edge
			int x_i = edges_[getEdgeNodeIndex(iEdge, 0)];
			int x_j = edges_[getEdgeNodeIndex(iEdge, 1)];

			// Add pairwise terms
			if (curWeight0 > 1e-8 || curWeight1 > 1e-8)
				graphCut -> add_edge(x_i, x_j, curWeight0, curWeight1);
		}

		// Run the GC solver
		tic = clock();
		double currentEnergy = graphCut -> maxflow();
		flowTime += (clock() - tic) / CLOCKS_PER_SEC;

		// Update the energy and the solution if the current configuration is good
		if (currentEnergy + hSum < energy_) {
			energy_ = currentEnergy + hSum;

			// Get the solution
			labeling_.resize(nodeNumber_, -1);
			for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
				labeling_[iNode] = graphCut -> what_segment(iNode);
			}
		}
		if (verbosityLevel_ >= 2) 
			printf("Iteration: %d; energy: %f\n", iIteration, currentEnergy + hSum + energyConstant_);

		if (iIteration < orderOfClassSwiches.size()) { // if not the last iteration prepare for the next
			// Make the switch to the next iteration
			classSwitches[orderOfClassSwiches[iIteration]] = 1 - classSwitches[orderOfClassSwiches[iIteration]];

			// update the hSum
			if(classSwitches[orderOfClassSwiches[iIteration]] == 1)
				hSum -=  classTreshold_[orderOfClassSwiches[iIteration]] * (1 - alpha_);
			else
				hSum +=  classTreshold_[orderOfClassSwiches[iIteration]] * (1 - alpha_);

			// cleat the graph 
			graphCut -> reset();
		}

	}

	if (verbosityLevel_ >= 1) 
		printf("Flow time: %fs\n", flowTime);

	// Get the lower bound
	lowerBound_ = energy_;

	// Set the result flag	
	resultsComputed_ = true;
	numNodesUnlabeled_ = 0;

	delete graphCut;
}


void CoopCutMinimizer::createHanoiTowerArray(int N, vector<int> *answer) const
{
	if (N == 1) {
		answer -> resize(1);
		answer -> at(0) = 0;
	} else {
		createHanoiTowerArray(N - 1, answer);
		int curSize = (int) answer -> size();
		answer -> push_back(N - 1);
		for(int i = curSize - 1; i >= 0; --i)
			answer -> push_back(answer -> at(i));
	}
}
void CoopCutMinimizer::setOrderOfClassSwitches(const vector<double>& classPriority, vector<int>* switches) const
{
	//construct and array - of priority and index
	vector<PairDoubleInt> pairArray(classPriority.size());
	for(int iClass = 0; iClass < classPriority.size(); ++iClass) {
		pairArray[iClass].first = classPriority[iClass];
		pairArray[iClass].second = iClass;
	}
	//sort the array of pairs
	std::sort(pairArray.begin(), pairArray.end(), less_PairDoubleInt);

	// find how many elements are with negative priority
	int lowPrCounter = 0;
	while (lowPrCounter < classPriority.size() && pairArray[lowPrCounter].first < 0) ++lowPrCounter;

	// Check the case when nothing should be checked
	if(lowPrCounter == classPriority.size()) {
		switches -> resize(0);
		return;
	}

	//Create a Hanoi tower array
	createHanoiTowerArray((int)classPriority.size() - lowPrCounter, switches);

	// update answer according to the initial sorting
	for(int iPos = 0; iPos < switches -> size(); ++iPos) {
		switches -> at(iPos) = pairArray[switches -> at(iPos) + lowPrCounter].second;
	}
}



void CoopCutMinimizer::minimize_exhaustive_GraphCut_dynamic()
{
	//Determine priority of classes for an exhaustive search: -1.0 - not included, largest number - the most expensive change
	vector<double> classPriority(edgeClassNumber_, -1.0);
	for (int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		if (classTreshold_[iClass] < (1.0 - (1e-2)) * totalEdgeClassWeight_[iClass]) { // class is not modular
			classPriority[iClass] = totalEdgeClassWeight_[iClass];
		}

	// create a switch order for an exhaustive search
	vector<int> orderOfClassSwiches(0, 0);
	setOrderOfClassSwitches(classPriority, &orderOfClassSwiches);

	// inialize the starting switch configuration: each switch corresponds to h_g
	vector<int> classSwitches(edgeClassNumber_, 1); // Class that are modular are always on

	// compute the sum of current values of h_g
	double hSum = 0;
	for(int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		hSum += (1 - classSwitches[iClass]) * classTreshold_[iClass] * (1 - alpha_);

	//create a graphCut object
	GraphCutObject* graphCut = new GraphCutObject(nodeNumber_, edgeNumber_);

	// make the first cut
	// Add nodes
	graphCut -> add_node(nodeNumber_);
  
	// add unary potentials for x_i
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		graphCut -> add_tweights(iNode, unaryPotentials_[iNode], 0);
	}

	// Get the arcs of each class
	vector<vector<GraphCutObject::arc_id>> arcsPerClass(edgeClassNumber_, vector<GraphCutObject::arc_id>(0, NULL));
	vector<vector<double>> arcAlternatingWeightsPerClass(edgeClassNumber_, vector<double>(0, 0));
	int node0, node1;
	GraphCutObject::arc_id arcPtr = graphCut -> get_first_arc();
	
	//compute time only for maxFlow
	double flowTime = 0;
	double tic = 0;


	//Add edges
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		//find the classes of current edges
		int curClass0 = edgeClass_[getEdgeClassIndex(iEdge, 0)];
		int curClass1 = edgeClass_[getEdgeClassIndex(iEdge, 1)];

		//find the weight of the current edge in the current direction
		double curWeight0 = (classSwitches[curClass0] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 0)] : edgeWeights_[getEdgeWeightIndex(iEdge, 0)] * alpha_;
		double curWeight1 = (classSwitches[curClass1] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 1)] : edgeWeights_[getEdgeWeightIndex(iEdge, 1)] * alpha_;
		
		// find the indices of all nodes related to the current edge
		int x_i = edges_[getEdgeNodeIndex(iEdge, 0)];
		int x_j = edges_[getEdgeNodeIndex(iEdge, 1)];

		// Add pairwise terms
		graphCut -> add_edge(x_i, x_j, curWeight0, curWeight1);

		// save references to both arcs
		arcsPerClass[curClass0].push_back(arcPtr);
		arcAlternatingWeightsPerClass[curClass0].push_back(edgeWeights_[getEdgeWeightIndex(iEdge, 0)] * (1 - alpha_));
		arcPtr = graphCut -> get_next_arc(arcPtr);
			
		arcsPerClass[curClass1].push_back(arcPtr);
		arcAlternatingWeightsPerClass[curClass1].push_back(edgeWeights_[getEdgeWeightIndex(iEdge, 1)] * (1 - alpha_));
		arcPtr = graphCut -> get_next_arc(arcPtr);
	}

	// Run the GC solver
	tic = clock();
	energy_ = graphCut -> maxflow() + hSum;
	flowTime += (clock() - tic) / CLOCKS_PER_SEC;

	// Get the solution  of unlabeled nodes
	labeling_.resize(nodeNumber_, -1);
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		labeling_[iNode] = graphCut -> what_segment(iNode);
	}

	// flow constant for reparametrizations
	double overflow = 0;
		
	for(int iIteration = 0; iIteration < orderOfClassSwiches.size(); ++iIteration) {
	
		//detect the current switch
		int curSwitch = orderOfClassSwiches[iIteration];

		// Make the switch to the next iteration
		classSwitches[curSwitch] = 1 - classSwitches[curSwitch];

		// update the hSum
		hSum += classTreshold_[curSwitch] * (1.0 - alpha_) * ((classSwitches[curSwitch] == 1)? -1 : 1);
		
		// update the graph for the marginals
		overflow += updateTheGraph(graphCut, classSwitches[curSwitch], arcsPerClass[curSwitch], arcAlternatingWeightsPerClass[curSwitch]); 

		// run the maxflow
		tic = clock();
		double curFlow = graphCut -> maxflow(true) - overflow;
		flowTime += (clock() - tic) / CLOCKS_PER_SEC;

		if (curFlow + hSum < energy_) {
			energy_ = curFlow + hSum;

			// Get the solution
			labeling_.resize(nodeNumber_, -1);
			for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
				labeling_[iNode] = graphCut -> what_segment(iNode);
			}
		}
		
		if (verbosityLevel_ >= 2) 
			printf("Iteration: %d; energy: %f\n", iIteration, curFlow + hSum + energyConstant_);
	}

	if (verbosityLevel_ >= 1) 
		printf("Flow time: %fs\n", flowTime);

	// Get the lower bound
	lowerBound_ = energy_;

	// Set the result flag	
	resultsComputed_ = true;

	// Set the number of unlabeled nodes
	numNodesUnlabeled_ = 0;
	
	delete graphCut;
}


double CoopCutMinimizer::updateTheGraph(GraphCutObject* graphCut, int switchDirection, const vector<GraphCutObject::arc_id> &arcs, 	const vector<double> &weights) const
{
	double overflow = 0;
	// update all the edges of this type
	for(int iEdge = 0; iEdge < arcs.size(); ++iEdge) {
		GraphCutObject::arc_id curArc = arcs[iEdge];
		double curChange =	weights[iEdge];
		double curResidual = graphCut -> get_rcap(curArc);

		// find the touched nodes
		int node0, node1;
		graphCut -> get_arc_ends(curArc, node0, node1);
		graphCut -> mark_node(node0);
		graphCut -> mark_node(node1);
	
		if (switchDirection == 1){
			graphCut -> set_rcap(curArc, curResidual + curChange);
		} else {
		if (curResidual - curChange > 0)
			graphCut -> set_rcap(curArc, curResidual - curChange);
		else {
			// do reparametrization
			double backFlowVal = -curResidual + curChange;
			
			//set the current edge
			graphCut -> set_rcap(curArc, 0);

			// adjust the sister arc
			//double tmp = graphCut -> get_rcap(curArc -> sister);
			graphCut -> set_rcap(curArc -> sister, graphCut -> get_rcap(curArc -> sister)  - backFlowVal);

			// adjust unary nodes
			graphCut -> add_tweights(node0, backFlowVal, 0);
			graphCut -> add_tweights(node1, 0, backFlowVal);

			// adjust an overflow
			overflow += backFlowVal;
		}
		}
	}
	return overflow;
}


void CoopCutMinimizer::minimize_greedy_marginals(int maxIter)
{
	//Determine priority of classes for an exhaustive search: -1.0 - not included, largest number - the most expensive change
	vector<double> classPriority(edgeClassNumber_, -1.0);
	for (int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		if (classTreshold_[iClass] < (1.0 - (1e-2)) * totalEdgeClassWeight_[iClass]) { // class is not modular
			classPriority[iClass] = totalEdgeClassWeight_[iClass];
		}

	// create a switch order for an exhaustive search
	//vector<int> orderOfClassSwiches(0, 0);
	//setOrderOfClassSwitches(classPriority, &orderOfClassSwiches);

	// inialize the starting switch configuration: each switch corresponds to h_g
	vector<int> classSwitches(edgeClassNumber_, 1); // Class that are modular are always on

	// compute the sum of current values of h_g
	double hSum = 0;
	for(int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		hSum += (1 - classSwitches[iClass]) * classTreshold_[iClass]  * (1 - alpha_);

	//create a graphCut object
	GraphCutObject* graphCut = new GraphCutObject(nodeNumber_, edgeNumber_);

	// make the first cut
	// Add nodes
	graphCut -> add_node(nodeNumber_);
  
	// add unary potentials for x_i
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		graphCut -> add_tweights(iNode, unaryPotentials_[iNode], 0);
	}

	// Get the arcs of each class
	vector<vector<GraphCutObject::arc_id>> arcsPerClass(edgeClassNumber_, vector<GraphCutObject::arc_id>(0, NULL));
	vector<vector<double>> arcAlternatingWeightsPerClass(edgeClassNumber_, vector<double>(0, 0));
	int node0, node1;
	GraphCutObject::arc_id arcPtr = graphCut -> get_first_arc();
	
	//compute time only for maxFlow
	double flowTime = 0;
	double tic = 0;

	//Add edges
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		//find the classes of current edges
		int curClass0 = edgeClass_[getEdgeClassIndex(iEdge, 0)];
		int curClass1 = edgeClass_[getEdgeClassIndex(iEdge, 1)];

		//find the weight of the current edge in the current direction
		double curWeight0 = (classSwitches[curClass0] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 0)] : edgeWeights_[getEdgeWeightIndex(iEdge, 0)] * alpha_;
		double curWeight1 = (classSwitches[curClass1] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 1)] : edgeWeights_[getEdgeWeightIndex(iEdge, 1)] * alpha_;
		
		// find the indices of all nodes related to the current edge
		int x_i = edges_[getEdgeNodeIndex(iEdge, 0)];
		int x_j = edges_[getEdgeNodeIndex(iEdge, 1)];

		// Add pairwise terms
		graphCut -> add_edge(x_i, x_j, curWeight0, curWeight1);

		// save references to both arcs
		arcsPerClass[curClass0].push_back(arcPtr);
		arcAlternatingWeightsPerClass[curClass0].push_back(edgeWeights_[getEdgeWeightIndex(iEdge, 0)] * (1.0 - alpha_));
		arcPtr = graphCut -> get_next_arc(arcPtr);
			
		arcsPerClass[curClass1].push_back(arcPtr);
		arcAlternatingWeightsPerClass[curClass1].push_back(edgeWeights_[getEdgeWeightIndex(iEdge, 1)] * (1.0 - alpha_));
		arcPtr = graphCut -> get_next_arc(arcPtr);
	}

	// Run the GC solver
	tic = clock();
	energy_ = graphCut -> maxflow() + hSum;
	flowTime += (clock() - tic) / CLOCKS_PER_SEC;

	// Get the solution
	labeling_.resize(nodeNumber_, -1);
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		labeling_[iNode] = graphCut -> what_segment(iNode);
	}

	// flow constant for reparametrizations
	double overflow = 0;
	
	vector<double> alternativeMarginals(edgeClassNumber_, 0);
	
		
	for(int iIteration = 0; iIteration < maxIter; ++iIteration) {

		//compute alternative marginals for al classes
		for(int iClass = 0; iClass < edgeClassNumber_; ++iClass) {
			// check that the current class is not modular
			if (classPriority[iClass] < 0) { continue; }
			
			
			// update the hSum
			double newHSum = hSum + classTreshold_[iClass] * (1.0 - alpha_) * ((1 - classSwitches[iClass] == 1)? -1 : 1);
			
			// update the graph for the marginals
			overflow += updateTheGraph(graphCut, 1 - classSwitches[iClass], arcsPerClass[iClass], arcAlternatingWeightsPerClass[iClass]); 

			// run the maxflow
			tic = clock();
			alternativeMarginals[iClass] = graphCut -> maxflow(true) - overflow + newHSum;
			flowTime += (clock() - tic) / CLOCKS_PER_SEC;

			// set everything back
			// update all the edges of this type
			overflow += updateTheGraph(graphCut, classSwitches[iClass], arcsPerClass[iClass], arcAlternatingWeightsPerClass[iClass]); 
		}

		int bestSwitch = -1;
		double bestEnergyDrop = 0;
		for(int iClass = 0; iClass < edgeClassNumber_; ++iClass) {
			// check that the current class is not modular
			if (classPriority[iClass] < 0) { continue; }
			
			if( energy_ - alternativeMarginals[iClass] > bestEnergyDrop) {
				bestEnergyDrop = energy_ - alternativeMarginals[iClass];
				bestSwitch = iClass;
			}
		}

		if (bestSwitch == -1) {
			break;
		} else {
			energy_ = energy_ - bestEnergyDrop;
	
			// update the switch
			classSwitches[bestSwitch] = 1 - classSwitches[bestSwitch];

			// update all the edges of this type
			overflow += updateTheGraph(graphCut, classSwitches[bestSwitch], arcsPerClass[bestSwitch], arcAlternatingWeightsPerClass[bestSwitch]); 

			// update hSum
			hSum += classTreshold_[bestSwitch] * (1.0 - alpha_) * ((classSwitches[bestSwitch] == 1)? -1 : 1);
			
		}

		if (verbosityLevel_ >= 2)
			printf("Iteration: %d; energy: %f\n", iIteration, energy_ + energyConstant_);
	}

	// compute the final maxflow
	energy_ = graphCut -> maxflow(true) - overflow + hSum;

	// Get the solution  of unlabeled nodes
	labeling_.resize(nodeNumber_, -1);
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		labeling_[iNode] = graphCut -> what_segment(iNode);
	}

	if (verbosityLevel_ >= 1) 
		printf("Flow time: %fs\n", flowTime);

	// Get the lower bound
	lowerBound_ = energy_;

	// Set the result flag	
	resultsComputed_ = true;

	// Set the number of unlabeled nodes
	numNodesUnlabeled_ = 0;
	
	delete graphCut;
}

void CoopCutMinimizer::minimize_allsteps_marginals(int maxIter)
{
	//Determine priority of classes for an exhaustive search: -1.0 - not included, largest number - the most expensive change
	vector<double> classPriority(edgeClassNumber_, -1.0);
	for (int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		if (classTreshold_[iClass] < (1.0 - (1e-2)) * totalEdgeClassWeight_[iClass]) { // class is not modular
			classPriority[iClass] = totalEdgeClassWeight_[iClass];
		}

	// create a switch order for an exhaustive search
	//vector<int> orderOfClassSwiches(0, 0);
	//setOrderOfClassSwitches(classPriority, &orderOfClassSwiches);

	// inialize the starting switch configuration: each switch corresponds to h_g
	vector<int> classSwitches(edgeClassNumber_, 1); // Class that are modular are always on

	// compute the sum of current values of h_g
	double hSum = 0;
	for(int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		hSum += (1 - classSwitches[iClass]) * classTreshold_[iClass]  * (1 - alpha_);

	//create a graphCut object
	GraphCutObject* graphCut = new GraphCutObject(nodeNumber_, edgeNumber_);

	// make the first cut
	// Add nodes
	graphCut -> add_node(nodeNumber_);
  
	// add unary potentials for x_i
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		graphCut -> add_tweights(iNode, unaryPotentials_[iNode], 0);
	}

	// Get the arcs of each class
	vector<vector<GraphCutObject::arc_id>> arcsPerClass(edgeClassNumber_, vector<GraphCutObject::arc_id>(0, NULL));
	vector<vector<double>> arcAlternatingWeightsPerClass(edgeClassNumber_, vector<double>(0, 0));
	int node0, node1;
	GraphCutObject::arc_id arcPtr = graphCut -> get_first_arc();
	
	//compute time only for maxFlow
	double flowTime = 0;
	double tic = 0;

	//Add edges
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		//find the classes of current edges
		int curClass0 = edgeClass_[getEdgeClassIndex(iEdge, 0)];
		int curClass1 = edgeClass_[getEdgeClassIndex(iEdge, 1)];

		//find the weight of the current edge in the current direction
		double curWeight0 = (classSwitches[curClass0] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 0)] : edgeWeights_[getEdgeWeightIndex(iEdge, 0)] * alpha_;
		double curWeight1 = (classSwitches[curClass1] == 1)? edgeWeights_[getEdgeWeightIndex(iEdge, 1)] : edgeWeights_[getEdgeWeightIndex(iEdge, 1)] * alpha_;
		
		// find the indices of all nodes related to the current edge
		int x_i = edges_[getEdgeNodeIndex(iEdge, 0)];
		int x_j = edges_[getEdgeNodeIndex(iEdge, 1)];

		// Add pairwise terms
		graphCut -> add_edge(x_i, x_j, curWeight0, curWeight1);

		// save references to both arcs
		arcsPerClass[curClass0].push_back(arcPtr);
		arcAlternatingWeightsPerClass[curClass0].push_back(edgeWeights_[getEdgeWeightIndex(iEdge, 0)] * (1.0 - alpha_));
		arcPtr = graphCut -> get_next_arc(arcPtr);
			
		arcsPerClass[curClass1].push_back(arcPtr);
		arcAlternatingWeightsPerClass[curClass1].push_back(edgeWeights_[getEdgeWeightIndex(iEdge, 1)] * (1.0 - alpha_));
		arcPtr = graphCut -> get_next_arc(arcPtr);
	}

	// Run the GC solver
	tic = clock();
	energy_ = graphCut -> maxflow() + hSum;
	flowTime += (clock() - tic) / CLOCKS_PER_SEC;

	// flow constant for reparametrizations
	double overflow = 0;
	
	vector<double> alternativeMarginals(edgeClassNumber_, 0);
		
	for(int iIteration = 0; iIteration < maxIter; ++iIteration) {
		bool changed = false;
		//compute alternative marginals for al classes
		for(int iClass = 0; iClass < edgeClassNumber_; ++iClass) {
			// check that the current class is not modular
			if (classPriority[iClass] < 0) { continue; }
			
			
			// update the hSum
			double newHSum = hSum + classTreshold_[iClass] * (1.0 - alpha_) * ((1 - classSwitches[iClass] == 1)? -1 : 1);
			
			// update the graph for the marginals
			overflow += updateTheGraph(graphCut, 1 - classSwitches[iClass], arcsPerClass[iClass], arcAlternatingWeightsPerClass[iClass]); 

			// run the maxflow
			tic = clock();
			alternativeMarginals[iClass] = graphCut -> maxflow(true) - overflow + newHSum;
			flowTime += (clock() - tic) / CLOCKS_PER_SEC;

			if (alternativeMarginals[iClass] < energy_) {
				// update the energy
				energy_ = alternativeMarginals[iClass];
				//update the switch
				classSwitches[iClass] = 1 - classSwitches[iClass];
				// update hSum
				hSum = newHSum;
				
				changed = true;

				if (verbosityLevel_ >= 2) 
					printf("Iteration: %d; class: %d; energy: %f\n", iIteration, iClass, energy_ + energyConstant_);
			} else {
				// set everything back
				// update all the edges of this type
				overflow += updateTheGraph(graphCut, classSwitches[iClass], arcsPerClass[iClass], arcAlternatingWeightsPerClass[iClass]); 
			}
		}
		if (!changed) break;
	}

	// compute the final maxflow
	energy_ = graphCut -> maxflow(true) - overflow + hSum;

	// Get the solution  of unlabeled nodes
	labeling_.resize(nodeNumber_, -1);
	for (int iNode = 0; iNode < nodeNumber_; ++iNode) {
		labeling_[iNode] = graphCut -> what_segment(iNode);
	}
	// The current energy is not necessarily the minimum over h, thus resompute the energy
	energy_ = computeEnergy(labeling_);


	if (verbosityLevel_ >= 1) 
		printf("Flow time: %fs\n", flowTime);

	
	// Set the result flag	
	resultsComputed_ = true;

	// Set the number of unlabeled nodes
	numNodesUnlabeled_ = 0;
	
	delete graphCut;
}


double CoopCutMinimizer::computeEnergy(const vector<int>& labeling) const
{
	double energy = 0;

	// add the unaries
	for(int iNode = 0; iNode < nodeNumber_; ++iNode) {
		if (labeling[iNode] == 1)
			energy += unaryPotentials_[iNode];
	}

	// compute the cut weight for each class of edges
	vector<double> edgeClassCut(edgeClassNumber_, 0.0);
	for(int iEdge = 0; iEdge < edgeNumber_; ++iEdge) {
		// edge forward
		if ( labeling[edges_[getEdgeNodeIndex(iEdge, 0)]] == 0 && labeling[edges_[getEdgeNodeIndex(iEdge, 1)]] == 1)
			edgeClassCut[edgeClass_[getEdgeClassIndex(iEdge, 0)]] += edgeWeights_[getEdgeWeightIndex(iEdge, 0)];

		//edge backward
		if ( labeling[edges_[getEdgeNodeIndex(iEdge, 0)]] == 1 && labeling[edges_[getEdgeNodeIndex(iEdge, 1)]] == 0)
			edgeClassCut[edgeClass_[getEdgeClassIndex(iEdge, 1)]] += edgeWeights_[getEdgeWeightIndex(iEdge, 1)];
	}

	// sum up the cooperation potentials
	for(int iClass = 0; iClass < edgeClassNumber_; ++iClass)
		energy += (edgeClassCut[iClass] < classTreshold_[iClass]) ? edgeClassCut[iClass] : classTreshold_[iClass] + alpha_ * (edgeClassCut[iClass] - classTreshold_[iClass]);

	return energy;
}

