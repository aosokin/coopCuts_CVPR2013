function [correctAnswer, presentColors] = prepareData2( imageName)

%imageName = '1_11_s';

image = im2double( imread( fullfile( 'images', [imageName, '.bmp'] ) ));

groundtruth = im2double(imread( fullfile( 'labelings', [imageName, '.bmp'] ) ));

pixelColors = reshape( groundtruth, [size(image, 1) * size(image, 2), 3] );

presentColors = unique(pixelColors, 'rows');

nodeNumber = size(image, 1) * size(image, 2);
classNumber = size(presentColors, 1);

%% create unary file
% read the seed file
seedImage = imread(fullfile('seeds', [imageName, '.png']));
seedColors = unique(reshape(seedImage, [], 3), 'rows');
seeds = false(size(image, 1), size(image, 2), classNumber);
unaryTerms = zeros(size(image, 1), size(image, 2), classNumber);
correctAnswer = nan(size(image, 1), size(image, 2));

for iColor = 1 : size(seedColors, 1)
    if seedColors(iColor, 1) ~= 255 || seedColors(iColor, 2) ~= 255 || seedColors(iColor, 3) ~= 255
       
        curMask = seedImage(:, :, 1) == seedColors(iColor, 1) & seedImage(:, :, 2) == seedColors(iColor, 2) & seedImage(:, :, 3) == seedColors(iColor, 3);
        
        for iLabel = 1 : classNumber
            curClassMask = groundtruth(:, :, 1) == presentColors(iLabel, 1) & groundtruth(:, :, 2) == presentColors(iLabel, 2) & groundtruth(:, :, 3) == presentColors(iLabel, 3);
            correctAnswer(curClassMask) = iLabel;
            
            seeds(:, :, iLabel) = seeds(:, :, iLabel) | (curMask & curClassMask);
        end
    end
end
for iLabel = 1 : classNumber
    mask = seeds(:, :, iLabel);
    tmp = makeUnaryPotentials( image, mask ); 
    unaryTerms(:, :, iLabel) = tmp;
end
for iLabel = 1 : classNumber
    mask = seeds(:, :, iLabel) == 1;
    for jLabel = 1 : classNumber
        if iLabel ~= jLabel
            tmp = unaryTerms(:, :, jLabel);
            tmp(mask) = 1000;
            unaryTerms(:, :, jLabel) = tmp;
        end
    end
end

unaryData = reshape(unaryTerms, [nodeNumber, classNumber])';

unaryData = unaryData / 50;

if ~exist('binaryFiles', 'dir')
    mkdir('binaryFiles');
end

unaryFileName = fullfile('binaryFiles', [imageName, '_unary.bin']);
fp = fopen(unaryFileName,'w');
fwrite(fp, int32(nodeNumber), 'int32');
fwrite(fp, int32(classNumber), 'int32');
fwrite(fp, unaryData(:), 'double');
fclose(fp);

%% create pairwise file
edgeClassNumber = 10;

[ edges, edgeWeights, edgeClasses ] = makePairwisePotentials( image, edgeClassNumber );
% 0-based indexing; edgeClassNumber - not included in any class

pairwiseFileName = fullfile('binaryFiles', [imageName, '_pairwise.bin']);
fp = fopen( pairwiseFileName, 'w' );
fwrite( fp, int32( size(edges, 1) ), 'int32' );
fwrite( fp, int32( edgeClassNumber ), 'int32' );
for iEdge = 1 : size(edges, 1)
    fwrite( fp, int32( edges(iEdge, 1) - 1 ), 'int32' );
    fwrite( fp, int32( edges(iEdge, 2) - 1), 'int32' );
    fwrite( fp, double( edgeWeights( iEdge ) ), 'double' );
    fwrite( fp, int32( edgeClasses(iEdge, 1) ), 'int32' );
    fwrite( fp, int32( edgeClasses(iEdge, 2) ), 'int32' );
end
fclose(fp);




