#include <stdio.h>

#include "CoopCutMultilabel.h"

/*
const char *unaryFile = "1_11_s_unary.bin";
const char *pairwiseFile = "1_11_s_pairwise.bin";
const char *outputFile = "output.bin";
*/

#include <time.h>

struct Options_coopCut {
	//input files
	char* pairwiseFile;
	char* unaryFile;
	char* outputFile;

	// method parameters;
	double edgeClassThreshold;
	double lambda;
	double alpha;
	int maxIter;
	
	// type of solver to use
	int solverType;
	
};

int main(int argc, char** argv)
{
	if (argc < 8) {
		printf("Usage: executableName UnaryFileName PairwiseFileName outputLabelFile lambda theta alpha solverType [maxIter]\n");
		printf("lambda - relative weight of pairwise terms\n");
		printf("theta - threshold for the break point: \\vartheta in Stefanie's paper\n");
		printf("alpha - relative coefficient of linear function after the breakpoint (should be < 1)\n");
		printf("solverType - type of optimization to run:\n");
		printf("\t0 - alpha-expansion a modular energy\n");
		printf("\t1 - greedy optimization over h with making all step that improve the energy\n");
		printf("\tmaxIter - maximum number of sweeps over all h (default: 10), only for solverType = 1\n");
		return 1;
	}
	
	Options_coopCut options;
	options.edgeClassThreshold = atof(argv[5]);
	options.lambda = atof(argv[4]);
	options.alpha = atof(argv[6]);
	
	options.pairwiseFile = argv[2];
	options.unaryFile = argv[1];
	options.outputFile = argv[3];

	options.solverType = atoi(argv[7]);
	
	options.maxIter = 10;
	if (argc > 8) options.maxIter = atoi(argv[8]);
	
	CoopCutMultilabel* coopCut = new CoopCutMultilabel();

	coopCut -> setAlpha(options.alpha);
	coopCut -> setLambda(options.lambda);
	coopCut -> setTheta(options.edgeClassThreshold);
	coopCut -> setGreedyMaxIter(options.maxIter);

	bool exitFlag = false;
	double tStart = clock();
	
	printf("Reading data: "); 
	exitFlag = coopCut -> readUnaryTermsFromFile( options.unaryFile );
	if (!exitFlag)  printf("ERROR while reading unary file\n");
	
	tStart = clock();
	exitFlag = coopCut -> readPairwiseTermsFromFile( options.pairwiseFile );
	if (!exitFlag) printf("ERROR while reading pairwise file\n");
	printf("Time: %f\n", (clock() - tStart) / CLOCKS_PER_SEC);


	int solverType = options.solverType;

	// switch the solve name
	char* methodName = NULL;
	switch(solverType) {
	case 0:
		methodName = "alpha-expansion";
		break;
	case 1:
		methodName = "greedy";
		break;
	default:
		printf("ERROR: unknown method specified!\n");
	}
	
	double elapsedTime;
	if (methodName != NULL) {
		//Run solver
		 tStart = (double)clock();
		//if (options.verbosityLevel >= 1)
			printf("Running %s...\n", methodName);
		
		switch (solverType) {
		case 0:
			exitFlag = coopCut -> minimizeEnergy_alphaExpansion();
			break;
		case 1:
			exitFlag = coopCut -> mimimizeEnergy_greedy();
			break;
		
		}
	
		elapsedTime = (clock() - tStart) / CLOCKS_PER_SEC;
		//if (options.verbosityLevel >= 1)
		if (!exitFlag)  
			printf("ERROR\n");
		else
			printf("%s finishes. Time: %fs; energy: %f\n", methodName, elapsedTime, coopCut -> getEnergy());
	} else {
		printf("ERROR: unknown method specified!\n");
		return 1;
	}
	
	printf("Writing results to file\n");
	exitFlag = coopCut -> writeSolutionToFile(options.outputFile);
	if (!exitFlag)  printf("ERROR while writing the output file\n");



	printf("Freeing memory\n");
	tStart = clock();
	delete coopCut;
	printf("Time: %f\n", (clock() - tStart) / CLOCKS_PER_SEC);

	// getchar();

	return 0;
}