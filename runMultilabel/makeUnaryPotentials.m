function [unaryTerms] = makeUnaryPotentials(image, mask)

A = double(image);

Av = reshape(A, [size(A, 1) * size(A, 2), 3] );

n1 = size(A,1);
n2 = size(A,2);
n = n1*n2;

disp('creating unary potentials');

objPix = (mask(:) == 1);
backgPix = false(size(objPix));
backgPix(1 : 10) = true;

% from www.calinon.ch/sourcecodes.php
% (adjust then the paths in the file makeGMM.m)

%if ~exist(sprintf('%s.gz',sname),'file') || overwrite

% k-means unaries
[logProb, tmp] = makeGMM(A, objPix, backgPix);
unaryTerms = reshape(logProb, [n1, n2]);

end



