function [unaryFile, edgeWeightFile, edgeClassFile] = makeBinFiles(imageFile, userLabelFile, numClasses, ...
    outputDir, createUnaries, createEdgeClasses, createEdgeWeights, unaryGMMflag)
%
% INPUT:
%   imageFile - the initial image file
%   userLabelFile - the file with user supplied labels
%   numClasses - the number of clusters for edge clustering
%   outputDir - the filder to write output binary files (default: 'temp')
%   createUnaries - flag for creating unary potentials: unaryFile (default: true)
%   createEdgeClasses - flag for creating edge weights: edgeWeightFile (default: true)
%   createEdgeWeights - flag for creating edge classes: edgeClassFile (default: true)
%   unaryGMMflag - flag for using GMM as unary potentials (default: true). false stands for bin histograms
%
% OUTPUT:
%   unaryFile - file with unary potentials
%   edgeWeightFile  - file with edge weights
%   edgeClassFile - file with edge classes
%
% This code is example code to compute unary potentials (here via
% histograms; one can also do GMMs) and edge classes needed to do image
% segmentation with cooperative cuts. It saves a file with unary potentials
% and one with edge classes.
% If log ratios are used, this only computes edge classes and not unary
% potentials.
%
% This needs the Matlab statistics toolbox.
%
% Code by Anton Osokin for CVPR 2013 paper
% "A Principled Deep Random Field Model for Image Segmentation"
%
% Code by Stefanie Jegelka and Jeff Bilmes for  CVPR 2011 paper "Submodularity beyong submodular
% energies: coupling edges in graph cuts" was heavily used
% 
% This code requires GMM fitting code by Sylvain Calinon
% http://www.calinon.ch/download/GMM-GMR-v2.0.zip
%
% Please acknowledge the corresponding papers in resulting publications.

% parameter initialization
edgeClassFile = [];
edgeWeightFile = [];
unaryFile = [];
if nargin < 8
    unaryGMMflag = true;
end
if nargin < 7
    createEdgeWeights = true;
end
if nargin < 6
    createEdgeClasses = true;
end
if nargin < 5
    createUnaries = true;
end
if nargin < 4
    outputDir = 'temp';
end


% get the name of the image
[~, imageName] = fileparts(imageFile);

% read in image
A = im2double(imread(imageFile));
N = size(A, 1);
M = size(A, 2);

% make set of pixels of the image
Av = reshape(A, [], 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create edge classes

% We start by taking differences of adjacent pixels. The order here is the
% same that is used in the C++ code provided online, it is important that
% the order of the entries is retained.
% This is done in four blocks (directions).

nodeId = reshape(1 : N * M, [N, M]);

% vertical
% create elist-part and use that to compute gradients
tmp1 = nodeId(1 : end - 1, :);
tmp2 = nodeId(2 : end, :);
edgeNodeIndex = [ tmp1(:), tmp2(:) ];

% horizontal, to right and to left
tmp1 = nodeId(:, 2 : end);
tmp2 = nodeId(:, 1 : end - 1);
edgeNodeIndex = [ edgeNodeIndex; tmp1(:), tmp2(:) ];

% diagonal \
tmp1 = nodeId(1 : end - 1, 1 : end - 1);
tmp2 = nodeId(2 : end, 2 : end);
edgeNodeIndex = [ edgeNodeIndex; tmp1(:), tmp2(:) ];

% diagonal /
tmp1 = nodeId(2 : end, 1 : end - 1);
tmp2 = nodeId(1 : end - 1, 2 : end);
edgeNodeIndex = [ edgeNodeIndex; tmp1(:), tmp2(:) ];

% make feature vectors for each edges
edgeFeatures = Av(edgeNodeIndex(:, 1), :) - Av(edgeNodeIndex(:, 2), :);

clear edgeNodeIndex

% now we are done creating the vectors of pixel differences
% next we transform it to the standard (Gaussian) edge weights
edgeWeights = sum( edgeFeatures .^ 2, 2);
nullWeights = (edgeWeights < eps);

sigma = mean( edgeWeights );
edgeWeights = 0.05 + 0.95 * exp( -edgeWeights / (2 * sigma));

if createEdgeWeights
    fprintf('Creating edge weights.\n')
    edgeWeightFile = [outputDir, '/', imageName, '_wts.bin'];
    
    % created uni-directed weights
    % these must be doubled in C code
    
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    fp = fopen(edgeWeightFile,'wb');
    fwrite(fp, N, 'int32'); %height
    fwrite(fp, M, 'int32'); %width
    fwrite(fp, length(edgeWeights), 'int32');
    fwrite(fp, edgeWeights, 'double');
    fclose(fp);
    
end

if createEdgeClasses
    fprintf('Creating edge classes.\n')
    
    % cluster the difference vectors into num_classes classes; the ones where
    % the difference is zero will go in an extra class
    
    edgeWeights = reshape([edgeWeights'; edgeWeights'], [], 1);
    nullWeights = reshape([nullWeights'; nullWeights'], [], 1);
    numWeights = length(edgeWeights);
    totalWeights = sum(edgeWeights);
    
    
    edgeFeatures = [ edgeFeatures, -edgeFeatures]';
    edgeFeatures = reshape(edgeFeatures, 3, [])';
    
    
    classes = ones(numWeights, 1) * numClasses;
    
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    edgeClassFile = [outputDir, '/',imageName, '_cl', num2str(numClasses), '.bin'];
    
    
    % k-means gives vector if assignments
    fprintf('Running k-means. This might need some time.\n')
    restcl = kmeans(edgeFeatures((nullWeights == 0), :), numClasses, 'EmptyAction', 'singleton','Distance', 'sqEuclidean');
    classes(nullWeights == 0) = restcl - 1;
    
    %%%%%%%%% WRITE EDGE CLASSES INTO A FILE %%%%%%%%%%%%%%%%%
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    fp = fopen(edgeClassFile,'wb');
    fwrite(fp, numClasses, 'int32');
    fwrite(fp, classes, 'int32');
    fclose(fp);
    
    for i2 = min(classes) : max(classes)
        fprintf('class %d: %d edges\n', i2, sum(classes == i2));
    end
    % end
    
end

clear edgeFeatures;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create unary potentials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if createUnaries
    fprintf('Creating unary potentials.\n')
    
    % read in the label file, coded by red and blue pixels
    L = im2double(imread( userLabelFile )) ;
    if size(L, 1) ~= size(A, 1) || size(L, 2) ~= size(A, 2)
        error(['Sizes of image file ', imageFile, ' and label file ', userLabelFile, 'do not match!' ]);
    end
    objPix = (L(:,:,1) > 0.5) & (L(:,:,2) < 0.5) & (L(:,:,3) < 0.5);
    backgPix = (L(:,:,1) < 0.5) & (L(:,:,1) < 0.5) & (L(:,:,3) > 0.5);
    objPix = find(objPix(:));
    backgPix = find(backgPix(:));
    clear L;
    
    if ~unaryGMMflag
        
        unaryFile = [outputDir, '/', imageName, '_ubinU.bin'];
        
        binsz = 4;
        % the next function returns two weight vectors, one for the source and
        % one for the sink edges
        [wo,wb] = make3Dhist(binsz, Av * 255, objPix, backgPix);
        wsink = -wo;
        wsource = -wb;
        
        % the pixels that the used marked receive heavy weigths
        hardwt = (totalWeights + sum(wsink))/2;
        wsource(objPix) = hardwt;
        wsink(backgPix) = hardwt;
        
        wsource = reshape(reshape(wsource, N, M)',[],1);
        wsink = reshape(reshape(wsink, N, M)',[],1);
        
        if ~exist(outputDir,'dir')
            mkdir(outputDir);
        end
        fp = fopen(unaryFile,'wb');
        fwrite(fp, wsource, 'double');
        fwrite(fp, wsink, 'double');
        fclose(fp);
        save(sprintf('%s-mat.mat',unaryFile),'wsource','wsink');
        
    else %%%% requires the GMM-GMR-v2.0 code
        % from www.calinon.ch/sourcecodes.php
        % (adjust then the paths in the file makeGMM.m)
        
        unaryFile = [outputDir, '/', imageName, '_ugmmU.bin'];
        
        % k-means unaries
        [wsource,wsink] = makeGMM(A * 255, backgPix, objPix);
        hardwt = (sum(totalWeights) + sum(wsource) + sum(wsink))/3;
        wsource(objPix) = hardwt;
        wsink(backgPix) = hardwt;
        wsource = reshape(reshape(wsource, N, M)',[],1);
        wsink = reshape(reshape(wsink, N, M)',[],1);
        
        if ~exist(outputDir,'dir')
            mkdir(outputDir);
        end
        fp = fopen(unaryFile,'wb');
        fwrite(fp, wsource, 'double');
        fwrite(fp, wsink, 'double');
        fclose(fp);
        
    end
    
    clear wsource;
    clear wsink;
    
end

end


