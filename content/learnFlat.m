function [template bx by] = learnFlat(folder,sxBySy,ax)
% clear

load config1
ay = ax;
nsketch = 100;

% (bx, by) should be adaptive per category

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
Ctransform(nimage, norient, allfiltered, imageSizes(:,1), imageSizes(:,2), Upperbound);

sx = bx; sy = by;
% local mean on saturated Gabor responses
localHalfx = ax; localHalfy = ax;
localAveFilter = ones(localHalfx,localHalfy,'single')/single(localHalfx)/single(localHalfy);
S1LocalMeanOverOrientation = cell(nimage,1);
for i = 1:nimage
	S1LocalMean = cell(1,norient);
	for o = 1:norient
		S1LocalMean{o} = filter2(localAveFilter,allfiltered{i,o},'same');
	end
	S1LocalMeanOverOrientation{i} = zeros(sx,sy);
	for o = 1:norient
		S1LocalMeanOverOrientation{i} = S1LocalMeanOverOrientation{i} + 1/norient * S1LocalMean{o};
    end
    S1LocalMeanOverOrientation{i} = Upperbound - S1LocalMeanOverOrientation{i};
end

template = learnFlatLogLinear(S1LocalMeanOverOrientation,nsketch,ax);
% =================================================================
save(['rawmodel_flat_' folder '_size' num2str(sxBySy) '_ax' num2str(ax) '.mat'],'template', 'bx', 'by', 'ax');


