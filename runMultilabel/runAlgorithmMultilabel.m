function curResult = runAlgorithmMultilabel(curImageFileName, ...
    curUnaryPotentialFileName, ...
    curPairwisePotenatialFileName, ...
    lambda, theta, alpha, ...
    algoInfo)

currentFileCode = randomString(20);
curResultFileName = ['result_', algoInfo, '_', num2str(round(lambda * (1e+5))), '_', num2str(round(theta * (1e+5))), '_', num2str(round(alpha * (1e+5))), '_', currentFileCode, '.bin'];


% Usage: executableName UnaryFileName PairwiseFileName outputLabelFile lambda theta alpha solverType [maxIter]
%     lambda - relative weight of pairwise terms
%     theta - threshold for the break point: \\vartheta in Stefanie's paper
%     alpha - relative coefficient of linear function after the breakpoint (should be < 1)
%     solverType - type of optimization to run:
%         0 - alpha-expansion a modular energy
%         1 - greedy optimization over h with making all step that improve the energy
% 	maxIter - maximum number of sweeps over all h (default: 10), only for solverType = 1


algoLine = '';
if strcmpi(algoInfo, 'alpha-expansion')
    algoLine = '0';
end
if strcmpi(algoInfo, 'greedy')
    algoLine = '1';
end
if strcmpi(algoInfo, 'greedy 1 iter')
    algoLine = '1 1';
end

commandLine =['coopcutMultilabel', ' ', ...
    '"', curUnaryPotentialFileName, '"', ' ', ...
    '"', curPairwisePotenatialFileName, '"',' ', ...
    '"', curResultFileName, '"', ' ', ...
    num2str(lambda), ' ', ...
    num2str(theta), ' ', ...
    num2str(alpha), ' ', ...
    algoLine];


system(commandLine, '-echo');
curResult = struct;
[curResult.labels, curResult.energy, curResult.time, curResult.methodEnergy] = analyzeResultsMultilabel(curResultFileName, curUnaryPotentialFileName, curPairwisePotenatialFileName, lambda, theta, alpha);

image = imread(curImageFileName);

curResult.imageLabels = reshape(curResult.labels, [size(image, 1), size(image, 2)] );

% clean up the temp file
delete(curResultFileName);
end

