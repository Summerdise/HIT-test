function get_natural_stat_flat(ax)


load config1
binsize = .1;  % binsize for computing F
Upperbound = 6.; % saturation value

negFolder = 'rand';  % image input folder
Iname = dir([negFolder '/*.jpg']);
nimage = min(length(Iname),200); % #images
Iname = Iname(1:nimage);
Isize = zeros(nimage, 2);
I = cell(1,nimage);
for i = 1 : nimage
    tmpIm = imread([negFolder '/' Iname(i).name]);
    if size(tmpIm,3) == 3
        tmpIm = rgb2gray(tmpIm);
    end
    I{i} = double(tmpIm); 
    Isize(i, :) = size(I{i});
end
sx = min(Isize(:, 1)); sy = min(Isize(:, 2)); 
for i = 1 : nimage
    I{i} = I{i}(1:sx, 1:sy); 
end   % make the sizes of the testing images the same


% FILTERING IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allfilter = makefilter(scale, norient);  % generate Gabor filters 
h = (size(allfilter{1}, 1)-1)/2;  % half size of Gabor
allfiltered = applyfilterfftsame(I, allfilter);  % filter training images
mex Ctransform.c
Ctransform(nimage, norient, allfiltered, Isize(:,1), Isize(:,2), Upperbound);

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
end


sub = 10;
count = 1;
binnum = floor(Upperbound/binsize)+1;
histog = zeros(binnum,1);
for i = 1:nimage
	tmp = S1LocalMeanOverOrientation{i}(1+h+ax:sub:sx-ax-h,1+h+ax:sub:sy-ax-h);
	tmp = Upperbound - tmp;
	current_histog = hist( tmp(:), binnum );
    histog = histog + current_histog';
end
histog = histog/sum(histog)/binsize;

% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = (0:(binnum-1))*binsize;
M = 50;
e = zeros(M, 1); 
Z = zeros(M, 1);
for k=1:M
    lambda = (k-1.)/10.;
    p = exp(lambda*r).*(histog');
    Z(k) = sum(p*binsize); p = p/Z(k);
    %plot(r, p);
    e(k) = sum(r.*p*binsize);
end
lam = (0:(M-1))/10.; 
lz = log(Z); 

save('nlf_flat.mat', 'M',...
    'lam', 'e', 'lz', 'Upperbound',...
    'scale', 'norient','ax','histog');

