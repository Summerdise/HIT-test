function [template bx by] = learnColor(folder,sxBySy,ax,gain_thres)

load config1
ay = ax;
nF = 20;
saturation = 0.001;

% (bx, by) should be adaptive per category

load storedColorHI storedHI

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
        I{i} = double(im); 
    else
        disp( 'not color image!' );
    	template = []; bx = 0; by = 0;
    	return
    end
end

bx = floor(mean(imageSizes(:,1)));
by = floor(mean(imageSizes(:,2)));
for i = 1:nimage
    I{i} = imresize(I{i},[bx by]);
end
imageSizes(:,1) = bx; imageSizes(:,2) = by;

% compute the locally normalized color histograms
localAveFilter = ones(ax,ay,'single')/single(ax)/single(ay);
scalar = 1.0 / 250;
for i = 1:nimage
	R = I{i}(:,:,1);
	G = I{i}(:,:,2);
	B = I{i}(:,:,3);
	R = filter2(localAveFilter,R,'same');
	G = filter2(localAveFilter,G,'same');
	B = filter2(localAveFilter,B,'same');
	%tot = R + G + B;
	%R = R ./ tot; G = G ./ tot; B = B ./ tot;
	I{i}(:,:,1) = R * scalar;
	I{i}(:,:,2) = G * scalar;
	I{i}(:,:,3) = B * scalar;
end

% compute mean of histogram map
mI = zeros(bx,by,3);
for i = 1:nimage
	mI(:,:,1) = mI(:,:,1) + 1.0/nimage * I{i}(:,:,1);
	mI(:,:,2) = mI(:,:,2) + 1.0/nimage * I{i}(:,:,2);
	mI(:,:,3) = mI(:,:,3) + 1.0/nimage * I{i}(:,:,3);
end

% compute variance of historam map
vI = zeros(bx,by);
for i = 1:nimage
	res = I{i} - mI;
	vI = vI + max(sum(res.^2,3), saturation ) * 1.0/double(nimage);
end

% CONSTANTS

varMapSaved = vI;

% pursuit the least-variance positions
inhibited = zeros(bx,by);
selectedMeanHist = zeros(nF,3);
maxValue = max(vI(:));
vI([1:ax,end-ax+1:end],:) = maxValue+1;
vI(:,[1:ax,end-ax+1:end]) = maxValue+1;
minVar = zeros(nF,1);
Mx = zeros(nF,1); My = zeros(nF,1);

for k = 1:nF
    [minVar(k) ind] = min(vI(:));
    
    
    col = ceil(ind/bx); My(k) = col;
    row = ind - (col-1)*bx; Mx(k) = row;
    % calculate the mean gradient histogram
    selectedMeanHist(k,:) = mI(row,col,:);
    
    % re-estimate the normalizing constant
    responses = sum((storedHI - repmat(selectedMeanHist(k,:),[size(storedHI,1),1])).^2,2);
    responses = max(responses,saturation);
    dZ(k) = mean(exp( - responses / 2 / minVar(k) ) / (sqrt(2*pi*minVar(k)))); % additional normalizing constant
    
    if -1/2-log( sqrt(2*pi.*minVar(k)) )-log(dZ(k)) < gain_thres || minVar(k) > maxValue
        nF = k-1;
        break;
    end
    
    % do inhibition
    left = max(1,floor(col-ay*1.5));
    right = min(by,floor(col+ay*1.5));
    top = max(1,floor(row-ax*1.5));
    bottom = min(bx,floor(row+ax*1.5));
    inhibited(top:bottom,left:right) = 1;
    vI(top:bottom,left:right) = maxValue + 1;
end
minVar = minVar(1:nF); Mx = Mx(1:nF); My = My(1:nF); selectedMeanHist = selectedMeanHist(1:nF,:); 
gain = (-1/2-log(sqrt(2*pi.*minVar))) - log(dZ(1:nF)');

template = struct('var',minVar,'meanHist',selectedMeanHist,...
    'inhibited',inhibited,'varMap',varMapSaved,'Mx',Mx,'My',My...
    , 'ax',ax,'ay',ay,'sx',bx,'sy',by,'gain',gain,'dZ',dZ, 'mI', mI);

% =================================================================
save(['rawmodel_color_' folder '_size' num2str(sxBySy) '_ax' num2str(ax) '.mat'],'template', 'bx', 'by', 'ax');


