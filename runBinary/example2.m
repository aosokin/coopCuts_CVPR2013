% The example of using methods presented in
% P. Kohli, A. Osokin, S. Jegelka. "A Principled Deep Random Field Model for Image Segmentation", 
% IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2013. 
%
% This example computes unary potentials, edge weights and edge classes

% dataset from CVPR 2011: S. Jegelka and J. Bilmes. "Submodularity beyond submodular energies: coupling edges in graph cuts"
imageName = {'001_oct', '002_scarf', '003_comb', '004_key', '005_bee', '006_insectjaw', ...
                '007_insectjaw', '008_mouse', '009_lamp', '010_blackbeetle', '012_greybeetle', ...
                '014_bamboo', '017_darktree', '019_redtree1', '020_smallplant', '028_darklamp' };
                
% The image
imageIndex = 5;
curImageFileName = fullfile('dataset', 'images', [imageName{imageIndex}, '.jpg']);
curGroundtruthFile = fullfile('dataset', 'groundtruth', [imageName{imageIndex}, '.png']);
curSeedFile = fullfile('dataset', 'seeds', [imageName{imageIndex}, '.png']);

% prepare data
outputDir = 'temp';
numClasses = 10;
[ curUnaryPotentialFileName, curEdgeWeightsFileName, curEdgeClassFileName ] = makeBinFiles(curImageFileName, curSeedFile, numClasses, outputDir);

% run algorithms
lambda = 1.5;
theta = 50 / 10 ^ 4;
alpha = 0.00;

algoInfo = {'GC_modular', ... % Regular graph cut that ignores high-order terms
            'coopcode', ... % Stefanie's method (CVPR 2011)
            'GC_global_dynamic', ... % Global minimum with dynamic cuts and search order heuristics
            'GC_greedy', ... % Greedy algorithm
            'GC_allSteps', ... % Greedy algorithm that makes a step if it improves the energy
            'GC_allSteps_1iter'}; % 1 iteration of the previous method
        
results = struct('labels', {}, 'energy', {}, 'time', {}, 'methodEnergy', {}, 'resultImage', {});        
for iAlgo = 1 : length(algoInfo)
    results(iAlgo) = runAlgorithm(curImageFileName, ...
                                    curUnaryPotentialFileName, ...
                                    curEdgeClassFileName, ...
                                    curEdgeWeightsFileName, ...
                                    lambda, theta, alpha, ...
                                    algoInfo{iAlgo});
                                
    if abs(results(iAlgo).methodEnergy -  results(iAlgo).energy) > (1e-5) * abs(results(iAlgo).energy)
        warning(['Method ', algoInfo{iAlgo}, ': something wrong with energy computation!']);
    end
    
    figure(iAlgo);
    imshow(results(iAlgo).resultImage);
end


