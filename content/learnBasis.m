function [template bx by] = learnBasis(folder,sxBySy)
% clear
load config1
nsketch = 80;
% folder = 'clockLHI';

% =================================================================
% Input training images
% -----------------------------------------------------------------

imageNames = dir([folder '/*.jpg']);
if isempty(imageNames)
    imageNames = dir([folder '/*.bmp']);
end
if isempty(imageNames)
    imageNames = dir([folder '/*.png']);
end
nimage = length(imageNames);
I = cell(1,nimage);
imageSizes = zeros(nimage,2);
for i = 1 : nimage
    im = imread([folder '/' imageNames(i).name]);
    sx = floor(sqrt(size(im,1)/size(im,2)*sxBySy));
    sy = floor(sx/size(im,1)*size(im,2));
    imageSizes(i,1) = sx; imageSizes(i,2) = sy;
    if size(im,3)==3
        im = rgb2gray(im);
    end
    I{i} = double(im);
end
% (bx, by) should be adaptive per category
bx = floor(mean(imageSizes(:,1)));
by = floor(mean(imageSizes(:,2)));
for i = 1:nimage
    I{i} = imresize(I{i},[bx by]);
end
imageSizes(:,1) = bx; imageSizes(:,2) = by;

% =================================================================
% FILTERING IMAGES
% -----------------------------------------------------------------
allfiltered = applyfilter(I, allfilter);  % filter training images
% -----------------------------------------------------------------
% WHITENING
% -----------------------------------------------------------------
Upperbound = 6;
Ctransform(nimage, norient, allfiltered, imageSizes(:,1), imageSizes(:,2), Upperbound);

dataWeight = ones(1,nimage)/nimage;
template = ABlearn2(dataWeight,allfiltered,allsymbol,allfilter,T,Lrange,Orange,sub,nsketch);

% =================================================================
save(['rawmodel_basis_' folder '_size' num2str(sxBySy) '.mat'],'template', 'bx', 'by');


