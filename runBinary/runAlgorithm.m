function curResult = runAlgorithm(curImageFileName, ...
    curUnaryPotentialFileName, ...
    curEdgeClassFileName, ...
    curEdgeWeightsFileName, ...
    lambda, theta, alpha, ...
    algoInfo)

currentFileCode = randomString(20);
curResultFileName = ['result_', algoInfo, '_', num2str(round(lambda * (1e+5))), '_', num2str(round(theta * (1e+5))), '_', num2str(round(alpha * (1e+5))), '_', currentFileCode, '.bin'];


if strcmpi(algoInfo, 'coopcode')
    % coopcode edgeClassFile edgeWeightsFile unaryFile lambda theta [outputLabelFile] [alpha]\n
    commandLine =['coopcode', ' ', ...
        curEdgeClassFileName, ' ', ...
        curEdgeWeightsFileName, ' ', ...
        curUnaryPotentialFileName, ' ', ...
        num2str(lambda), ' ', ...
        num2str(theta), ' ', ...
        curResultFileName, ' ', ...
        num2str(alpha)];
else
    % Usage: executableName edgeClassFileName edgeWeightsFileName UnaryFileName outputLabelFileName lambda theta alpha solverType [verbosityLevel] [maxIter]
    % lambda - relative weight of pairwise terms
    % theta - threshold for the break point: \vartheta in Stefanie's paper
    % alpha - relative coefficient of linear function after the breakpoint (should be < 1)
    % solverType - type of optimization to run:
    % 	0 - graph cut on a modular energy
    % 	1 - QPBO on a full energy with h and z
    % 	2 - exhaustive search over h without reusing the flow
    % 	3 - exhaustive search over h with reusing the flow
    % 	4 - greedy optimization over h
    % 	5 - greedy optimization over h with making all step that improve the energy
    % verbosityLevel - 0 - no printing; 1 - some printing (default); 2 - iteration printing
    % maxIter - maximum number of sweeps over all h (default: 10), only for solverType = 4, 5
    algoLine = '';
    if strcmpi(algoInfo, 'GC_modular')
        algoLine = '0 1';
    end
    if strcmpi(algoInfo, 'QPBO')
        algoLine = '1 1';
    end
    if strcmpi(algoInfo, 'GC_global')
        algoLine = '2 1';
    end
    if strcmpi(algoInfo, 'GC_global_dynamic')
        algoLine = '3 1';
    end
    if strcmpi(algoInfo, 'GC_greedy')
        algoLine = '4 1 10';
    end
    if strcmpi(algoInfo, 'GC_allSteps')
        algoLine = '5 1 10';
    end
    if strcmpi(algoInfo, 'GC_allSteps_1iter')
        algoLine = '5 1 1';
    end
    
    commandLine =['coopCutMinimizer', ' ', ...
        curEdgeClassFileName, ' ', ...
        curEdgeWeightsFileName, ' ', ...
        curUnaryPotentialFileName, ' ', ...
        curResultFileName, ' ', ...
        num2str(lambda), ' ', ...
        num2str(theta), ' ', ...
        num2str(alpha), ' ', ...
        algoLine];
    
end

[status, result] = system(commandLine, '-echo');

if status == 0
    % everything worked fine
    if exist(curResultFileName, 'file')
        % result file exists
        curResult = struct;
        [curResult.labels, curResult.energy, curResult.time, curResult.methodEnergy] = analyzeResults(curResultFileName, curUnaryPotentialFileName, curEdgeWeightsFileName, curEdgeClassFileName, lambda, theta, alpha);
        curResult.resultImage = makeResultImage(curImageFileName, curResult.labels);
        % clean up the temp file
        delete(curResultFileName);
    else
        error(['File "', curResultFileName, '" was not found. Check that command line "', commandLine, '" creates this file']);
    end
else
    if exist(curResultFileName, 'file')
        delete(curResultFileName);
    end
    error(['Command line "', commandLine, '" exited with nonzero exit code. Try running this line form a command line."']);
end

end

function string = randomString(stringLength)

symbols = ['a':'z', '0':'9'];
nums = randi(numel(symbols),[1 stringLength]);
string = symbols (nums);

end

function resultImage = makeResultImage(imageFileName, labels)

I = im2double(imread(imageFileName));
height = size(I, 1);
width = size(I, 2);

if size(labels, 1) ~=  height || size(labels, 2) ~=  width
    error('Image and segmentation are not compatible!');
end

resutlImage = nan(height, width, 3);
for iChannel = 1 : 3
    tmp = I(:, :, iChannel);
    
    tmp( labels == 0) = 1;
    tmp(labels < 0) = 0.5;
    
    resultImage(:, :, iChannel) = tmp;
end

end


