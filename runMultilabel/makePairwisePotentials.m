function [edges, edgeWeights, edgeClasses] = makePairwisePotentials(image, edgeClassNumber)

% read in image
A = double(image);

n1 = size(A,1);
n2 = size(A,2);
n = n1*n2;

Av = reshape(double(A), [], 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create edge classes


% We start by taking differences of adjacent pixels. The order here is the
% same that is used in the C++ code provided online, it is important that
% the order of the entries is retained.
% This is done in four blocks (directions).

disp('creating edges ...');

% node indices
edges = zeros(0, 2);
nodeIndices = reshape(1 : n, [n1, n2]);

% vertical
% create elist-part and use that to compute gradients
tmpli = [];
elist = [];
for i=1:n2
    tmpli = [tmpli; ...
        [([1:(n1-1)]'+(i-1)*n1), ([2:n1]'+(i-1)*n1)] ];
end
wts = Av(tmpli(:,1),:) - Av(tmpli(:,2),:);
edges = [edges; reshape(nodeIndices(1 : end - 1, :), [], 1), reshape(nodeIndices(2 : end , :), [], 1) ];
elist = [elist; tmpli];


% horizontal, from right and to left
tmpli = [];
for i=1:(n2-1)
    tmpli = [tmpli; ...
        [([1:n1]'+i*n1), ([1:n1]'+(i-1)*n1)] ];
end
wts1 = Av(tmpli(:,1),:) - Av(tmpli(:,2),:);
wts = [wts; [wts1]];
edges = [edges; reshape(nodeIndices(:, 2 : end), [], 1), reshape(nodeIndices(:, 1 : end - 1), [], 1) ];

elist = [elist; tmpli];

diagstart = size(elist,1)+1;
% diagonal \
tmpli = [];
for i=1:(n2-1)
    tmpli = [tmpli; ...
        [([1:(n1-1)]'+(i-1)*n1), ([2:n1]'+(i)*n1)] ];
end
wts1 = Av(tmpli(:,1),:) - Av(tmpli(:,2),:);
wts = [wts; [wts1]];
edges = [edges; reshape(nodeIndices(1 : end - 1, 1 : end - 1), [], 1), reshape(nodeIndices(2 : end, 2 : end), [], 1) ];
elist = [elist; tmpli];


% diagonal /
tmpli = [];
for i=1:(n2-1)
    tmpli = [tmpli; ...
        [([2:n1]'+(i-1)*n1), ([1:(n1-1)]'+i*n1)] ];
end
wts1 = Av(tmpli(:,1),:) - Av(tmpli(:,2),:);
wts = [wts; wts1];
edges = [edges; reshape(nodeIndices(2 : end, 1 : end - 1), [], 1), reshape(nodeIndices(1 : end - 1, 2 : end), [], 1) ];
elist = [elist; tmpli];



% now we are done creating the vectors of pixel differences
% next we transform it to the standard (Gaussian) edge weights
tmpli = sum( wts.^2, 2);
tmpw = reshape([tmpli'; tmpli'], [], 1);
nullinds = (tmpw < eps);

fprintf('\n%d zero differences\n', sum(nullinds));

sigma = mean(tmpli);
wtsfull = wts;
weights = 0.05+0.95*exp( - tmpli / (2*sigma));
edgeWeights = weights;


mw = max(weights);
m = size(elist,1);
clear tmpw;


% cluster the difference vectors into num_classes classes; the ones where
% the difference is zero will go in an extra class

disp('clustering edges ...');

% double each entry in wts, in negative form to match tmpw
wts = [wtsfull, -wtsfull]';
wts = reshape(wts, 3, [])';
clear wtsfull; clear elist;

classes = ones(2 * m,1) * edgeClassNumber;

distance = 'sqEuclidean';

% make edge classes by clustering (k-means)
str1 = RandStream.create('mrg32k3a', 'Seed', 27, 'NumStreams',1);
RandStream.setGlobalStream(str1);
% if ~exist(sprintf('%s.gz',sname),'file')

% k-means gives vector if assignments
restcl = kmeans( wts((nullinds==0),:), edgeClassNumber, 'EmptyAction', 'singleton', 'Distance', distance );
classes(nullinds==0) = restcl-1;

edgeClasses = reshape(classes, 2, [])';

end

