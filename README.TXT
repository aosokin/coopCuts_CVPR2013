This software library implements the algorithms for cooperative cuts described in
P. Kohli, A. Osokin, S. Jegelka. "A Principled Deep Random Field Model for Image Segmentation", 
IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2013. 

Primary contact:  Anton Osokin (anton.osokin@gmail.com)

All the implemented algorithms use the min-cut/max-flow algorithm code described in: 
http://pub.ist.ac.at/~vnk/software/maxflow-v3.02.src.tar.gz
Y. Boykov and V. Kolmogorov. "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision". 
In IEEE Transactions on Pattern Analysis and Machine Intelligence, September 2004.

The code heavily relies on the dynamic cut method:
P. Kohli and P.H.S. Torr. "Efficiently Solving Dynamic Markov Random Fields Using Graph Cuts."
International Conference on Computer Vision (ICCV), 2005 

If you use this software for research purposes, you should cite
the aforementioned papers in any resulting publication.

This package also includes
1) QPBO code:
http://pub.ist.ac.at/~vnk/software/QPBO-v1.3.src.tar.gz
V. Kolmogorov and C. Rother. "Minimizing non-submodular functions with graph cuts - a review". IEEE TPAMI, 29(7):1274-1279, July 2007.
2) The initial algorithm for the cooperative cuts:
http://melodi.ee.washington.edu/~jegelka/cc/code/coopcode.zip
S. Jegelka and J. Bilmes. "Submodularity beyond submodular energies: coupling 
edges in graph cuts". CVPR, 2011. 
3) Code for fitting Gaussian mixture models
http://www.calinon.ch/download/GMM-GMR-v2.0.zip
Calinon, S. (2009) Robot Programming by Demonstration: A Probabilistic Approach. EPFL/CRC Press.
Calinon, S., Guenter, F. and Billard, A. (2007) On Learning, Representing and Generalizing a Task in a Humanoid Robot. IEEE Transactions on Systems, Man and Cybernetics, Part B. 37:2, 286-298.

LICENSE & DISCLAIMER
---------------------

    Copyright 2013 Pushmeet Kohli (pkohli@microsoft.com), Anton Osokin (anton.osokin@gmail.com), Stefanie Jegelka (stefje@eecs.berkeley.edu)

    This software can be used for research purposes for free. In any
    resulting publication, you should cite papers mentioned above. For
    commercial use, please contact the authors above. 

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

USING THE CODE
-----------------------------

The code is supplied as MSVC 2010 solution (tested on win7, x64) with MATLAB scripts that call the main routines:

1) ./srcBinary/coopCutMinimizer - all the proposed optimization methods for binary energies
./srcBinary/coopCutMinimizer: main.cpp, coopCutMinimizer.h, coopCutMinimizer.cpp, minimizer.cpp

main.cpp provides a command line interface:
./coopCutMinimizer edgeClassFile edgeWeightsFile UnaryFile outputLabelFile lambda theta alpha solverType [verbosityLevel] [maxIter]
edgeClassFile, edgeWeightsFile, UnaryFile - input files with data (format by S. Jegelka)
outputLabelFile - file to save the results
lambda - relative weight of pairwise terms
theta - threshold for the break point: \\vartheta in Stefanie's paper
alpha - relative coefficient of linear function after the breakpoint (should be < 1)
solverType - type of optimization to run:
	0 - graph cut on a modular energy
	1 - QPBO on a full energy with h and z
	2 - exhaustive search over h without reusing the flow
	3 - exhaustive search over h with reusing the flow
	4 - greedy optimization over h
	5 - greedy optimization over h with making all step that improve the energy
verbosityLevel - 0 - no printing; 1 - some printing (default); 2 - iteration printing
maxIter - maximum number of sweeps over all h (default: 10), only for solverType = 4, 5

./srcBinary/graphCuts - the max-flow library by V. Kolmogorov (v.3.02)
src/graphCuts/maxflow-v3.02.src: graph.h, block.h, graph.cpp, maxflow.cpp, instances.inc

./srcBinary/QPBO - the QPBO library by Vladimir Kolmogorov (v.1.3)
src/QPBO/QPBO-v1.3.src:  QPBO.h, block.h, QPBO.cpp, QPBO_extra.cpp, QPBO_maxflow.cpp, QPBO_postprocessing.cpp, instances.inc

./srcBinary/coopcode - the initial coopcut algorithm by Stefanie Jegelka
src/coopcode: call.cpp, costfun.h, edge_data.h, edat_read.h, itBM.h, costfun.cpp, costfun_fulltrunc.cpp, edat_read.cpp, itBM.cpp

2) ./runBinary - Matlab scripts for running the code for binary energies
./runBinary/example1.m - the example RUNME file. Uses precomputed binary files from ./runBinary/data_example1
./runBinary/example2.m - the example RUNME file. Constructs binary files form raw images from ./runBinary/dataset
./runBinary/analyzeResults.m - recomputes the energy of the produced solution
./runBinary/runAlgorithm.m - calls c++ code from a console
./runBinary/makeBinFiles.m - constructs binary files with unary potentials, edge weights, and edge classes
./runBinary/coopCutMinimizer.exe - an executable for coopCutMinimizer (built with MSVC 2010, x64)
./runBinary/coopcode.exe - an executable for coopcode (built with MSVC 2010, x64)
./runBinary/data_example1 - data for exampleBee_preparedData.m
./runBinary/dataset - dataset form the paper
./runBinary/GMM-GMR-v2.0 - code by Sylvain Calinon for fitting GMMs

3) ./srcMultilabel/coopcutMultilabel - all the proposed optimization methods for multi-label energies
./srcMultilabel/coopCutMinimizer: main.cpp, CoopCutMultilabel.h, CoopCutMultilabel.cpp

To compile the code, please download the alpha-expansion code http://vision.csd.uwo.ca/code/gco-v3.0.zip and extract it into ./srcMultilabel/gcoLibrary/

main.cpp provides a command line interface:
./coopcutMultilabel UnaryFileName PairwiseFileName outputLabelFile lambda theta alpha solverType [maxIter]
UnaryFileName, PairwiseFileName - input binary files containing unary potentials and edge class data
outputLabelFile - name of a file  to save results (generated by this code)
lambda - relative weight of pairwise terms
theta - threshold for the break point: \\vartheta in Stefanie's paper
alpha - relative coefficient of linear function after the breakpoint (should be < 1)
solverType - type of optimization to run:
	0 - alpha-expansion a modular energy
	1 - greedy optimization over h with making all step that improve the energy
maxIter - maximum number of sweeps over all h (default: 10), works only for solverType = 1

4) ./runMultilabel - Matlab scripts for running the code for multi-label energies
./runMultilabel/exampleMultilabel.m - the example RUNME file. Construct binary files from images and runs the algorithm.
./runMultilabel/results - results generated for our paper


Note that the supplied binaries are build with MSVC 2010 and require MS .NET Framework 4:
http://www.microsoft.com/en-us/download/details.aspx?id=17851

