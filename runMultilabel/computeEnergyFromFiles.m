function energy = computeEnergyFromFiles( labels, unaryFile, pairwiseFile, lambda, theta, alpha)

% read the unary potentials
fid = fopen(unaryFile, 'r');
nodeNumber = fread(fid, 1, 'int32');
labelNumber = fread(fid, 1, 'int32');
unaryPotentials = fread(fid, labelNumber * nodeNumber, 'double');
fclose(fid);

unaryPotenatials = reshape( unaryPotentials, [labelNumber, nodeNumber] );

% read the pairwise potentials
fid = fopen( pairwiseFile, 'r' );
edgeNumber = fread( fid, 1, 'int32' );
edgeClassNumber = fread( fid, 1, 'int32' );

edges = nan(edgeNumber, 2);
edgeWeights = nan(edgeNumber, 1);
edgeClasses = nan(edgeNumber, 2);

for iEdge = 1 : edgeNumber
    edges(iEdge, 1) = fread( fid, 1, 'int32' ) + 1;
    edges(iEdge, 2) = fread( fid, 1, 'int32' ) + 1;
    edgeWeights(iEdge) = fread( fid, 1, 'double' );
    edgeClasses(iEdge, 1) = fread( fid, 1, 'int32' );
    edgeClasses(iEdge, 2) = fread( fid, 1, 'int32' );
end
fclose(fid);

% % read the label file
% fid = fopen( solutionFile, 'r' );
% curNodeNumber = fread(fid, 1, 'int32');
% % labels = double(fread(fid, curNodeNumber, 'schar')) + 1;
% time = fread(fid, 1, 'double');
% methodEnergy = fread(fid, 1, 'double');
% fclose(fid);

% compute the total weights of all classes
totalSum = zeros(edgeClassNumber, 1);
for iEdgeClass = 1 : edgeClassNumber
    mask = edgeClasses(:, 1) == (iEdgeClass - 1);
    totalSum(iEdgeClass) = sum(edgeWeights(mask));
    
    mask = edgeClasses(:, 2) == (iEdgeClass - 1);
    totalSum(iEdgeClass) = totalSum(iEdgeClass) + sum(edgeWeights(mask));
end

% compute the energy
unary = sum(unaryPotentials(labels + (0 : nodeNumber - 1)' * labelNumber));

pairwise = 0;

for iLabel = 1 : labelNumber
    for iEdgeClass = 0 : edgeClassNumber - 1
        % forward
        curEdges = (labels(edges(:, 1)) == iLabel) & (labels(edges(:, 2)) ~= iLabel) & (edgeClasses(:, 1) == iEdgeClass);
        curEdgeSum = 0.5 * sum(edgeWeights(curEdges));
        
        % backward
        curEdges = (labels(edges(:, 1)) ~= iLabel) & (labels(edges(:, 2)) == iLabel) & (edgeClasses(:, 2) == iEdgeClass);
        curEdgeSum = curEdgeSum + 0.5 * sum(edgeWeights(curEdges));
        b = (1 - alpha) * theta * totalSum(iEdgeClass + 1);
        pairwise = pairwise + min(curEdgeSum, alpha * curEdgeSum + b);

        
    end
end

pairwise = pairwise + sum((labels(edges(:, 1)) ~= labels(edges(:, 2))) .* (0.5 * edgeWeights) .* (edgeClasses(:, 1) ==  edgeClassNumber));
pairwise = pairwise + sum((labels(edges(:, 1)) ~= labels(edges(:, 2))) .* (0.5 * edgeWeights) .* (edgeClasses(:, 2) ==  edgeClassNumber));

energy = unary + lambda * pairwise;


end