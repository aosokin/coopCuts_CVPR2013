#include <stdio.h>
#include <vector>
using std::vector;

#include "coopCutMinimizer.h"
#include <time.h>

struct Options_coopCut {
	//input files
	char* edgeWeightFile;
	char* edgeClassFile;
	char* unaryFile;
	char* outputFile;

	// method parameters;
	double edgeClassThreshold;
	double edgeClassLastThreshold;
	double lambda;
	bool useDiagonalTerms;
	double alpha;
	int maxIter;
	int verbosityLevel;

	// type of solver to use
	int solverType;
	
};


// hierarchical cooperative cuts - example code using QPBO
int main(int argc, char** argv)
{
	if (argc < 9) {
		printf("Usage: executableName classFileName edgeWeightsFile UnaryFileName outputLabelFile lambda theta alpha solverType [verbosityLevel] [maxIter]\n");
		printf("lambda - relative weight of pairwise terms\n");
		printf("theta - threshold for the break point: \\vartheta in Stefanie's paper\n");
		printf("alpha - relative coefficient of linear function after the breakpoint (should be < 1)\n");
		printf("solverType - type of optimization to run:\n");
		printf("\t0 - graph cut on a modular energy\n");
		printf("\t1 - QPBO on a full energy with h and z\n");
		printf("\t2 - exhaustive search over h without reusing the flow\n");
		printf("\t3 - exhaustive search over h with reusing the flow\n");
		printf("\t4 - greedy optimization over h\n");
		printf("\t5 - greedy optimization over h with making all step that improve the energy\n");
		printf("verbosityLevel - 0 - no printing; 1 - some printing (default); 2 - iteration printing\n");
		printf("maxIter - maximum number of sweeps over all h (default: 10), only for solverType = 4, 5\n");
		return 1;
	}
	
	Options_coopCut options;
	options.edgeClassLastThreshold = 1.0;
	options.edgeClassThreshold = atof(argv[6]);
	options.lambda = atof(argv[5]);
	options.alpha = atof(argv[7]);
	options.useDiagonalTerms = true;

	options.edgeWeightFile = argv[2];
	options.edgeClassFile = argv[1];
	options.unaryFile = argv[3];
	options.outputFile = argv[4];

	options.solverType = atoi(argv[8]);
	
	options.verbosityLevel = 1;
	if (argc > 9) options.verbosityLevel = atoi(argv[9]);
	
	options.maxIter = 10;
	if (argc > 10) options.maxIter = atoi(argv[10]);
	
	//Options_coopCut options;
	//options.edgeClassLastThreshold = 1.0;
	//options.edgeClassThreshold = 0.003;
	//options.lambda = 1.0;
	//options.alpha = 0.1;
	//options.useDiagonalTerms = true;


	//options.edgeWeightFile = "testdata/smbee2_wts.bin";
	//options.edgeWeightFile = "testdata/beeNew_wts.bin";
	//options.edgeClassFile = "testdata/smbee2cl10.bin";
	//options.edgeClassFile = "testdata/beeNew_cl10nd.bin";
	//options.unaryFile = "testdata/smbee2ugmmU.bin";
	//options.unaryFile = "testdata/smbee2ubinU.bin";
	//options.unaryFile = "testdata/beeNew_ugmmU.bin";
	//options.outputFile = "testdata/smbee2_myseg.bin";
	//options.outputFile = "testdata/beeNew_myseg.bin";

	//options.solverType = 2;

  /*
  if (argc < 3){
    printf("use as: ./example_ho imname outfile-name \n, e.g. example_ho testdata/smbee testout/smbee_segment.bin");
    return 1;
  }
  */
  

  // input filename, e.g. "testdata/smbee"
  //const char* imname = "testdata/beeNew_";
  //const char* imname = "testdata/smbee";
  //const char* imname = argv[1];
  //const char* outfile = "testdata/beeNew_myseg.bin"; 
  //const char* outfile = "testdata/smbee_myseg.bin"; 
  //const char* outfile = argv[2];

	CoopCutMinimizer* solver = new CoopCutMinimizer();

	// set thresholds
	solver -> setClassTheta(options.edgeClassThreshold, options.edgeClassLastThreshold );

	// set lambda parameter
	solver -> setLambda(options.lambda);

	// set lambda parameter
	solver -> setAlpha(options.alpha);

	// set verbosity parameter
	solver -> setVerbosity(options.verbosityLevel);

	// read pairwise data
	if (solver -> readPairwisePotentials_Stefanie(options.edgeWeightFile, options.edgeClassFile, options.useDiagonalTerms)) {
		if (options.verbosityLevel >= 1)
			printf("Pairwise terms read successfully\n");
	} else {
		printf("Error while reading pairwise terms!\n");
		return 1;
	}

	// read unary data
	if (solver -> readUnaryPotentials_Stefanie(options.unaryFile)) {
		if (options.verbosityLevel >= 1)
			printf("Unary terms read successfully\n");
	} else {
		printf("Error while reading unary terms!\n");
		return 1;
	}

	int solverType = options.solverType;

	// switch the solve name
	char* methodName = NULL;
	switch(solverType) {
	case 0:
		methodName = "GC_modular";
		break;
	case 1:
		methodName = "QPBO";
		break;
	case 2: 
		methodName = "GC_exhaustive";
		break;
	case 3: 
		methodName = "GC_exhaustive_dynamic";
		break;
	case 4: 
		methodName = "GC_greedy";
		break;
	case 5: 
		methodName = "GC_greedy_allsteps";
		break;
	default:
		printf("ERROR: unknown method specified!\n");
	}
	double tStart;
	double elapsedTime;
	if (methodName != NULL) {
		//Run solver
		 tStart = (double)clock();
		if (options.verbosityLevel >= 1)
			printf("Running %s...\n", methodName);
		
		switch (solverType) {
		case 0:
			solver -> minimize_modular_GraphCut();
			break;
		case 1:
			solver -> minimize_QPBO();
			break;
		case 2:
			solver -> minimize_exhaustive_GraphCut();
			break;
		case 3:
			solver -> minimize_exhaustive_GraphCut_dynamic();
			break;
		case 4 :
			solver -> minimize_greedy_marginals(options.maxIter);
			break;
		case 5: 
			solver -> minimize_allsteps_marginals(options.maxIter);
			break;
		}
	
		elapsedTime = (clock() - tStart) / CLOCKS_PER_SEC;
		if (options.verbosityLevel >= 1)
			printf("%s finishes. Time: %fs; energy: %f; lower bound: %f; unlabeled: %d\n", methodName, elapsedTime, solver -> getEnergy(), solver -> getLowerBound(), solver -> getNumberNodesUnlabeled());
	} else {
		printf("ERROR: unknown method specified!\n");
		return 1;
	}


	//Write the labeling to file
	if (solver -> writeLabeling_Stefanie(options.outputFile, elapsedTime)) {
		if (options.verbosityLevel >= 1)
			printf("Labeling saved successfully\n");
	} else {
		printf("Error while saving the labeling!\n");
		return 1;
	}


	delete solver;

	//printf("Push any key\n");
	//getchar();

	return 0;
}
