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

