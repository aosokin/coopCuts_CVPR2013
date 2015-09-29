function [labels, energy, time, methodEnergy] = analyzeResults(labelFile, unaryFile, edgeWeightFile, edgeClassFile, lambda, theta, alpha)
% CAUTION: this version uses a counter-intuitive way to compute the cut: the edge is cut iff it goes from the object to the background

magicConstant = 50;

%% read the data
% read the edge weights
fid = fopen(edgeWeightFile, 'r');
height = fread(fid, 1, 'int32');
width = fread(fid, 1, 'int32');
edgeWeightNumber = fread(fid, 1, 'int32');
edgeWeights = fread(fid, edgeWeightNumber, 'double');
fclose(fid);

% read the edge classes
fid = fopen(edgeClassFile, 'r');
edgeClassNumber = fread(fid, 1, 'int32');
edgeClasses = fread(fid, edgeWeightNumber * 2, 'int32');
fclose(fid);

% read the unary potentials
fid = fopen(unaryFile, 'r');
unaryPotentials = fread(fid, 2 * width * height, 'double');
fclose(fid);

% read labels
fid = fopen(labelFile,'r');
N = fread(fid, 1, 'int32');
labels = fread(fid, N, 'schar');
time = fread(fid, 1, 'double');
methodEnergy = fread(fid, 1, 'double');
fclose(fid);

if height * width ~= N
    error('Image and segmentation are not compatible!');
end
labels = reshape(labels, [width, height])';

%% process the data to the trivial form
% the unaries
unaryTerms = nan(height, width, 2);
unaryTerms(:, :, 1) = reshape(unaryPotentials(1 : N), [width, height])';
unaryTerms(:, :, 2) = reshape(unaryPotentials(N + 1 : 2 * N), [width, height])';

% use the magic constant
edgeWeights = edgeWeights * magicConstant;

% the edge weights
offset = 0;
verticalEdgeWeights = reshape(edgeWeights((1 : (height - 1) * width) + offset), [height - 1, width]);
offset = offset + (height - 1) * width;
horizontalEdgeWeights = reshape(edgeWeights((1 : height * (width - 1)) + offset), [height, width - 1]);
offset = offset +  height * (width - 1);
diagonalMainEdgeWeights = reshape(edgeWeights((1 : (height - 1) * (width - 1)) + offset), [height - 1, width - 1]);
offset = offset +  (height - 1) * (width - 1);
diagonalSecondEdgeWeights = reshape(edgeWeights((1 : (height - 1) * (width - 1)) + offset), [height - 1, width - 1]);

% the edge classes
offset = 0;
verticalEdgeClasses = reshape(edgeClasses((1 : 2 * (height - 1) * width) + offset), [2, height - 1, width]);
verticalEdgeClasses = permute(verticalEdgeClasses, [2, 3, 1]);
offset = offset + 2 * (height - 1) * width;
horizontalEdgeClasses = reshape(edgeClasses((1 : 2 * height * (width - 1)) + offset), [2, height, width - 1]);
horizontalEdgeClasses = permute(horizontalEdgeClasses, [2, 3, 1]);
offset = offset +  2 * height * (width - 1);
diagonalMainEdgeClasses = reshape(edgeClasses((1 : 2 * (height - 1) * (width - 1)) + offset), [2, height - 1, width - 1]);
diagonalMainEdgeClasses = permute(diagonalMainEdgeClasses, [2, 3, 1]);
offset = offset +  2 * (height - 1) * (width - 1);
diagonalSecondEdgeClasses = reshape(edgeClasses((1 : 2 * (height - 1) * (width - 1)) + offset), [2, height - 1, width - 1]);
diagonalSecondEdgeClasses = permute(diagonalSecondEdgeClasses, [2, 3, 1]);

%% compute the energy
% compute the unary term
unaryTerm = sum(sum(unaryTerms(:, :, 1) .* (labels == 0), 2), 1);
unaryTerm = unaryTerm + sum(sum(unaryTerms(:, :, 2) .* (labels == 1), 2), 1);

% compute sum of edges cut and total sum of edges for all classes
edgeSumTotal = zeros(edgeClassNumber + 1, 1);
edgeSumCut = zeros(edgeClassNumber + 1, 1);
edgeNumberTotal = zeros(edgeClassNumber + 1, 1);
edgeNumberCut = zeros(edgeClassNumber + 1, 1);

for iClass = 1 : edgeClassNumber + 1
    % add vertical edges forward
    mask = verticalEdgeClasses(:, :, 1) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* verticalEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(1 : end - 1, :) == 1 & labels(2 : end, :) == 0;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* verticalEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));
    
    % add vertical edges backward
    mask = verticalEdgeClasses(:, :, 2) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* verticalEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(1 : end - 1, :) == 0 & labels(2 : end, :) == 1;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* verticalEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));

    % add horizontal edges forward   CAUTION: strange ordering here
    mask = horizontalEdgeClasses(:, :, 1) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* horizontalEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(:, 1 : end - 1) == 0 & labels(:, 2 : end) == 1;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* horizontalEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));

    % add horizontal edges backward
    mask = horizontalEdgeClasses(:, :, 2) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* horizontalEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(:, 1 : end - 1) == 1 & labels(:, 2 : end) == 0;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* horizontalEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));

    % add main diagonal edges '\' forward
    mask = diagonalMainEdgeClasses(:, :, 1) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* diagonalMainEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(1 : end - 1, 1 : end - 1) == 1 & labels(2 : end, 2 : end) == 0;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* diagonalMainEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));

    % add main diagonal edges '\' backward
    mask = diagonalMainEdgeClasses(:, :, 2) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* diagonalMainEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(1 : end - 1, 1 : end - 1) == 0 & labels(2 : end, 2 : end) == 1;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* diagonalMainEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));

    % add second diagonal edges '/' forward
    mask = diagonalSecondEdgeClasses(:, :, 1) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* diagonalSecondEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(2 : end, 1 : end - 1) == 1 & labels(1 : end - 1, 2 : end) == 0;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* diagonalSecondEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));

    % add second diagonal edges '/' backward
    mask = diagonalSecondEdgeClasses(:, :, 2) == iClass - 1;
    edgeSumTotal(iClass) = edgeSumTotal(iClass) + sum(mask(:) .* diagonalSecondEdgeWeights(:));
    edgeNumberTotal(iClass) = edgeNumberTotal(iClass) + sum(mask(:));
    maskCut = labels(2 : end, 1 : end - 1) == 0 & labels(1 : end - 1, 2 : end) == 1;
    edgeSumCut(iClass) = edgeSumCut(iClass) + sum(maskCut(:) .* mask(:) .* diagonalSecondEdgeWeights(:));
    edgeNumberCut(iClass) = edgeNumberCut(iClass) + sum(maskCut(:) .* mask(:));
end

% combine terms to compute the energy
energy = unaryTerm / lambda + edgeSumCut(edgeClassNumber + 1);
for iClass = 1 : edgeClassNumber
    b =  (1 - alpha) * theta * edgeSumTotal(iClass);
    energy = energy + min(edgeSumCut(iClass), alpha * edgeSumCut(iClass) + b);
end


end
