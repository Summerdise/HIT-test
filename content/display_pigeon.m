% function displayTemplateAndMatching

% Display active basis and adaboost templates;
% and display matched templates on example images.
clear
close all
nDisplay = 10;
folder = 'pigeonHead';
sxBySy = 14400;

% Change the following parameter to control how many features to show for the active basis template. 
sk_thres = 1.35;
tex_thres = 0.7;

%% active basis template from 5 training positives
load(sprintf('rawmodel_basis_%s_size%d',folder,sxBySy));
sx = template.sx; sy = template.sy;

load nlf
load config1
epsilon = .1;  % allowed correlation between selected Gabors
C = corr(allfilter, epsilon);    % generate the inhibition maps


sk_nF = sum( template.gain > sk_thres );
Mx = template.Mx(1:sk_nF);
My = template.My(1:sk_nF);
Mi = template.Mi(1:sk_nF);
Mm = template.Mm(1:sk_nF);
Mm1 = template.Mm1(1:sk_nF);

%% active basis template matching
Iname = dir( [folder '/*.jpg'] );
nimage = length(Iname); % # images
Isize = zeros(nimage, 2);
I = cell(1,nimage);
for i = 1 : nimage
    tmpIm = imread(sprintf('%s/%s',folder,Iname(i).name));
    if size(tmpIm,3) == 3
        tmpIm = rgb2gray(tmpIm);
    end
    I{i} = imresize(double(tmpIm),[sx sy]); 
    Isize(i, :) = size(I{i});
end

Sx = zeros(1, nimage) + sx; 
Sy = zeros(1, nimage) + sy; 

% filter the testing images 
disp('start filtering');
tic
allfiltered = applyfilterfftsame(I, allfilter);  % filter testing images
disp(['filtering time: ' num2str(toc) ' seconds']);

% sigmoid transformation
Ctransform(nimage, norient, allfiltered, Sx, Sy, Upperbound);

% prepare the spaces for C code
Asym = cell(1, nimage);   % symbolic plot for each image with active Gabor
for i = 1 : nimage
    Asym{i} = zeros(sx,sy);  
end
tss = zeros(nimage, sk_nF);  % total scores for testing images

% pass on the pointers to C for heavy computation
bkp = zeros(nimage,norient);
disp('start testing on sketch template');
tic
mex Ctesting.c
Ctesting(nimage, norient, allfiltered, C, h, sx, sy,...
    Mi, Mx, My, Mm, Mm1, allsymbol(1, :), Asym, Lrange,...
    Orange, sub, sk_nF, Upperbound, tss, M, lam, e, lz, bkp);


%%

K = nimage;
[tsssort, ind] = sort(tss);
J = cell(1,nimage*2+2);
J{1} = -sym;
J{2} = [];
for i = 1 : K 
    J{2*i+1} = I{ind(nimage-i+1)};
    J{2*i+2} = -(Asym{ind(nimage-i+1)});
end
towrite_sketch = displayImages(J(3:3+nDisplay-1),nDisplay,sx,sy);

% texture template matching
sxBySyTex = 150*150; ax = 10;
load(sprintf('rawmodel_tex_%s_size%d_ax%d',folder,sxBySyTex,ax));
tex_nF = sum(template.gain>tex_thres);
sx = template.sx; sy = template.sy;

I = cell(1,nimage);
for i = 1 : nimage
    tmpIm = imread(sprintf('%s/%s',folder,Iname(i).name));
    if size(tmpIm,3) == 3
        tmpIm = rgb2gray(tmpIm);
    end
    I{i} = imresize(double(tmpIm),[sx sy]); 
    Isize(i, :) = size(I{i});
end

Sx = zeros(1, nimage) + sx; 
Sy = zeros(1, nimage) + sy; 

% filter the testing images 
disp('start filtering');
tic
allfiltered = applyfilterfftsame(I, allfilter);  % filter testing images
disp(['filtering time: ' num2str(toc) ' seconds']);

% sigmoid transformation
Ctransform(nimage, norient, allfiltered, Sx, Sy, Upperbound);

% compute orientation histograms
hI = cell(size(allfiltered,1),size(allfiltered,2));
for i = 1:size(hI,1)
    for j = 1:size(hI,2)
        hI{i,j} = zeros(size(allfiltered{1,1}));
    end
end

% compute local mean maps
varMap = 1.5*ones(size(allfiltered{1,1}));
sub = max(floor(ax/5),floor(h/4));
CGraHist(nimage,norient,allfiltered,hI,varMap,Upperbound,sx,sy,ax,ax,h,sub);

disp('start testing texture template');
tic
nimage = size(hI,1);
norient = size(hI,2);

tss = zeros(nimage,tex_nF);
texAsym = cell(nimage,1);
selectedMeanHist = template.meanHist;
minVar = template.var;
dZ = template.dZ;
Mx = template.Mx;
My = template.My;
for i = 1:nimage
    texAsym{i} = zeros(sx,sy);
    for iF = 1:tex_nF
        tmpHI = zeros(1,norient);
        for o = 1:norient
            tmpHI(o) = hI{i,o}(Mx(iF),My(iF));
        end
        tss(i,iF) = -sum( (tmpHI-selectedMeanHist(iF,:)).^2) / 2 / (minVar(iF)+1e-10) - log(1e-10+sqrt(2*pi*minVar(iF))) - log(dZ(iF));
        if tss(i,iF) < 1
        	tss(i,iF) = 0;
        end
        vec = tmpHI;
        vec = vec.^4 + 1e-10;
        vec = vec/sum(vec);
        
        small_symbol = zeros( h*2+1, h*2+1 );
        for ori = 1:norient
        	small_symbol = small_symbol + vec(ori) * allsymbol{ori} * (tss(i,iF));
        end
        
        paved_symbol = zeros( size(small_symbol)*2 + 1 );
        paved_symbol(1:2*h+1,1:2*h+1) = small_symbol;
        paved_symbol(2*h+1+(1:2*h+1),1:2*h+1) = small_symbol;
        paved_symbol(1:2*h+1,2*h+1+(1:2*h+1)) = small_symbol;
        paved_symbol(2*h+1+(1:2*h+1),2*h+1+(1:2*h+1)) = small_symbol;
        paved_symbol = imresize(paved_symbol,[ax+2*h+1 ax+2*h+1]);
        row_ind = ( Mx(iF) - floor(size(paved_symbol,1)/2) -1 ) + (1:size(paved_symbol,1));
        col_ind = ( My(iF) - floor(size(paved_symbol,2)/2) -1 ) + (1:size(paved_symbol,2));
        texAsym{i}(row_ind,col_ind) = ...
            max(texAsym{i}(row_ind,col_ind) ...
            , paved_symbol );

    end
end

for i = 1 : nimage
    mixedSym = zeros(sx,sy,3);
    basisSym = -J{2*i+2};
    basisSym = imresize( basisSym, [sx, sy], 'bilinear' );
    mixedSym(:,:,2) = mixedSym(:,:,2) - basisSym ;
    mixedSym(:,:,1) = mixedSym(:,:,1) - basisSym ;
    mixedSym(:,:,3) = mixedSym(:,:,3) - basisSym ;
    texSym = texAsym{ind(nimage-i+1)};
    mixedSym(:,:,2) = mixedSym(:,:,2) - 1*texSym;
    mixedSym(:,:,3) = mixedSym(:,:,3) - 1*texSym;
    mixedSym = 255 * (mixedSym-min(mixedSym(:)))/(max(mixedSym(:))-min(mixedSym(:)));
    J{2*i+2} = mixedSym; 
end

towrite = displayImages(J(3:3+nDisplay-1),nDisplay,sx,sy);
imwrite(towrite,sprintf('matched_%s.png',folder));
% figure;imshow(towrite,[]);


