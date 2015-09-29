imageName = '1_11_s';

image = im2double( imread( fullfile( 'images', [imageName, '.bmp'] ) ));

groundtruth = im2double(imread( fullfile( 'labelings', [imageName, '.bmp'] ) ));

pixelColors = reshape( groundtruth, [size(image, 1) * size(image, 2), 3] );

presentColors = unique(pixelColors, 'rows');

for iColor = 1 : size(presentColors, 1)
   mask = groundtruth(:, :, 1) == presentColors(iColor, 1) & ...
            groundtruth(:, :, 2) == presentColors(iColor, 2) & ...
            groundtruth(:, :, 3) == presentColors(iColor, 3);
   mask = double( mask );
   
   imwrite(mask, fullfile( 'seeds', [imageName, '_class', num2str(iColor), '.png'] ));
            
end