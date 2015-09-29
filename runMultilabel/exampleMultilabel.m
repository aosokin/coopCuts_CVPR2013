% example

imageName = '1_11_s';

[correctAnswer, presentColors] = prepareData( imageName);

curImageFileName = ['images/',imageName, '.bmp'];
curUnaryPotentialFileName = ['binaryFiles/', imageName, '_unary.bin'];
curPairwisePotenatialFileName = ['binaryFiles/', imageName, '_pairwise.bin'];
lambda = 1.5;
theta = 50 / 10000;
alpha = 0;
%algoInfo = 'alpha-expansion';
%algoInfo = 'greedy';
algoInfo ='greedy 1 iter';

result = runAlgorithmMultilabel(curImageFileName, ...
    curUnaryPotentialFileName, ...
    curPairwisePotenatialFileName, ...
    lambda, theta, alpha, ...
    algoInfo);

imshow( result.imageLabels, [])