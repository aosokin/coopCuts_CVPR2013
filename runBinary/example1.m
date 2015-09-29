% The example of using methods presented in
% P. Kohli, A. Osokin, S. Jegelka. "A Principled Deep Random Field Model for Image Segmentation", 
% IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2013. 
%
% This example uses precomputed binary files stored in data_example1

% % Small image
% curImageFileName = 'data_example1/smbee2.jpg';
% curUnaryPotentialFileName = 'data_example1/smbee2ubinU.bin';
% curEdgeClassFileName = 'data_example1/smbee2cl10.bin';
% curEdgeWeightsFileName = 'data_example1/smbee2_wts.bin';

% % Real image
curImageFileName = 'data_example1/bee.jpg';
curUnaryPotentialFileName = 'data_example1/beeugmmU.bin';
curEdgeClassFileName = 'data_example1/beecl10.bin';
curEdgeWeightsFileName = 'data_example1/bee_wts.bin';

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


                                