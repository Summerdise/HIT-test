function tilt

% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = .7;  % scale of Gabors, length = 17 pixels
norient = 16;  % number of orientations
Upperbound = 6.; % saturation value
binsize = .2;  % binsize for computing F

% TRAINING EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
negFolder = 'rand';  % image input folder
negTrainImages = dir([negFolder '/*.jpg']);
nimage = length(negTrainImages); % #images
Isize = zeros(nimage, 2);
I = cell(1,nimage);
for i = 1 : nimage
    tmpIm = imread([negFolder '/' negTrainImages(i).name]);
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
disp('start filtering');
tic
allfilter = makefilter(scale, norient);  % generate Gabor filters 
h = (size(allfilter{1}, 1)-1)/2;  % half size of Gabor
allfiltered = applyfilterfftsame(I, allfilter);  % filter training images
disp(['filtering time: ' num2str(toc) ' seconds']);

binnum = floor(Upperbound/binsize)+1;  % binnumbers
histog = zeros(binnum, 1);  % store F  

% MEX-C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mex Chistog.c;   % compile C code
disp('start histogramming');
Chistog(nimage, norient, allfiltered, h, sx, sy, binsize, binnum, histog, Upperbound);
disp(['histogramming time: ' num2str(toc) ' seconds']);


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

save('nlf.mat', 'M',...
    'lam', 'e', 'lz', 'Upperbound',...
    'scale', 'norient');




