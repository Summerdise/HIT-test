function illustrateContrastTemplate

load config1
mex mexc_ComputeMAX1.cpp
pos = 'catHead';
posSx = 150; pos_thres = 1.2;
load(sprintf('rawmodel_basis_%s_size%d.mat',pos,posSx*posSx));
posT = template;
neg = 'bearHead'; negSx = 120; neg_thres = 0.8;
load(sprintf('rawmodel_basis_%s_size%d.mat',neg,negSx*negSx));
negT = template;
weight_thres = 0.7;
stepsize = 10;
[sym detailed_xyo] = contrastActiveBasis(posT,negT,allsymbol,weight_thres,pos_thres,neg_thres,stepsize);
imwrite(sym(:,:,1),sprintf('contrast_basis_%s_%s_pos.png',pos,neg));
imwrite(sym(:,:,2),sprintf('contrast_basis_%s_%s_com.png',pos,neg));
imwrite(sym(:,:,3),sprintf('contrast_basis_%s_%s_neg.png',pos,neg));
weight_thres = 0.01;
sym_svm = get_SVM_template(detailed_xyo,pos,neg,pos_thres,neg_thres,posSx,negSx,weight_thres,allfilter,allsymbol);
imwrite(sym_svm(:,:,1),sprintf('contrast_svm_%s_%s_pos.png',pos,neg));
imwrite(sym_svm(:,:,2),sprintf('contrast_svm_%s_%s_com.png',pos,neg));
imwrite(sym_svm(:,:,3),sprintf('contrast_svm_%s_%s_neg.png',pos,neg));
imwrite(1-posT.sym,sprintf('constrast_t%s.png',pos));
nSketchNeg = 50;
im = displayGaborTemplate([negSx,negSx], negT.Mx(1:nSketchNeg), negT.My(1:nSketchNeg), negT.Mi(1:nSketchNeg), ...
	ones(1,nSketchNeg), 0 );
im = 255 * (im-min(im(:)))/(max(im(:))-min(im(:)));
imwrite(im,sprintf('constrast_t%s.png',neg));
disp(neg)


pos = 'catHead';
sx = 150; pos_thres = 1.2;
load(sprintf('rawmodel_basis_%s_size%d.mat',pos,sx*sx));
posT = template;
neg = 'dogHead'; sx = 100; neg_thres = 1.2;
load(sprintf('rawmodel_basis_%s_size%d.mat',neg,sx*sx));
negT = template;
weight_thres = 0.8;
stepsize = 5;
[sym detailed_xyo] = contrastActiveBasis(posT,negT,allsymbol,weight_thres,pos_thres,neg_thres,stepsize);
imwrite(sym(:,:,1),sprintf('contrast_basis_%s_%s_pos.png',pos,neg));
imwrite(sym(:,:,2),sprintf('contrast_basis_%s_%s_com.png',pos,neg));
imwrite(sym(:,:,3),sprintf('contrast_basis_%s_%s_neg.png',pos,neg));
weight_thres = 0.01;
sym_svm = get_SVM_template(detailed_xyo,pos,neg,pos_thres,neg_thres,posSx,negSx,weight_thres,allfilter,allsymbol);
imwrite(sym_svm(:,:,1),sprintf('contrast_svm_%s_%s_pos.png',pos,neg));
imwrite(sym_svm(:,:,2),sprintf('contrast_svm_%s_%s_com.png',pos,neg));
imwrite(sym_svm(:,:,3),sprintf('contrast_svm_%s_%s_neg.png',pos,neg));
nSketchNeg = 50;
im = displayGaborTemplate([negSx,negSx], negT.Mx(1:nSketchNeg), negT.My(1:nSketchNeg), negT.Mi(1:nSketchNeg), ...
	ones(1,nSketchNeg), 0 );
im = 255 * (im-min(im(:)))/(max(im(:))-min(im(:)));
imwrite(im,sprintf('constrast_t%s.png',neg));
disp(neg)




pos = 'catHead';
sx = 150; pos_thres = 1.3;
load(sprintf('rawmodel_basis_%s_size%d.mat',pos,sx*sx));
posT = template;
neg = 'wolfHead'; sx = 120; neg_thres = 1.2;
load(sprintf('rawmodel_basis_%s_size%d.mat',neg,sx*sx));
negT = template;
weight_thres = 0.8;
stepsize = 5;
[sym detailed_xyo] = contrastActiveBasis(posT,negT,allsymbol,weight_thres,pos_thres,neg_thres,stepsize);
imwrite(sym(:,:,1),sprintf('contrast_basis_%s_%s_pos.png',pos,neg));
imwrite(sym(:,:,2),sprintf('contrast_basis_%s_%s_com.png',pos,neg));
imwrite(sym(:,:,3),sprintf('contrast_basis_%s_%s_neg.png',pos,neg));
weight_thres = 0.01;
sym_svm = get_SVM_template(detailed_xyo,pos,neg,pos_thres,neg_thres,posSx,negSx,weight_thres,allfilter,allsymbol);
imwrite(sym_svm(:,:,1),sprintf('contrast_svm_%s_%s_pos.png',pos,neg));
imwrite(sym_svm(:,:,2),sprintf('contrast_svm_%s_%s_com.png',pos,neg));
imwrite(sym_svm(:,:,3),sprintf('contrast_svm_%s_%s_neg.png',pos,neg));
imwrite(1-posT.sym,sprintf('constrast_t%s.png',pos));
nSketchNeg = 50;
im = displayGaborTemplate([negSx,negSx], negT.Mx(1:nSketchNeg), negT.My(1:nSketchNeg), negT.Mi(1:nSketchNeg), ...
	ones(1,nSketchNeg), 0 );
im = 255 * (im-min(im(:)))/(max(im(:))-min(im(:)));
imwrite(im,sprintf('constrast_t%s.png',neg));
disp(neg)

%%
function [basisSym detailed_xyo] = contrastActiveBasis(posT,negT,allsymbol,weight_thres,pos_thres, neg_thres, stepsize)
% weight_thres: the threshold for whether the feature weight is small
% for each scanned grid, only keep one basis, and accumulate the feature weight
h = floor( size( allsymbol{1} , 1 ) / 2 );
posScale = 120 / posT.sx;
negScale = 120 / negT.sx;
stepx = stepsize; stepy = stepsize; stepo = 2; norient = 16;
weights = zeros(120/stepx, (120/stepy), norient/2);
detailed_xyo = zeros(120/stepx, (120/stepy), norient/2,3);
% deal with positive template
nF = sum(posT.gain > pos_thres);
for iF = 1:nF
    ix = ceil(posT.Mx(iF)*posScale / stepx);
    iy = ceil(posT.My(iF)*posScale / stepy);
    io = ceil( (posT.Mi(iF)+1) / stepo );
    weights(ix,iy,io) = weights(ix,iy,io) + posT.Mm(iF);
    if detailed_xyo(ix,iy,io,1) == 0
        detailed_xyo(ix,iy,io,1:3) = [ceil(posT.Mx(iF)*posScale) ceil(posT.My(iF)*posScale) posT.Mi(iF)];
    end
end
% deal with negative template
nF = sum(negT.gain > neg_thres);
for iF = 1:nF
    found = false;
    for jF = 1:sum(posT.gain>pos_thres)
        if abs( posT.Mx(jF)*posScale - negT.Mx(iF)*negScale ) <= stepx &&...
                abs( posT.My(jF)*posScale - negT.My(iF)*negScale ) <= stepy &&...
                ( abs( posT.Mi(jF) - negT.Mi(iF) ) < stepo || abs( posT.Mi(jF) - negT.Mi(iF) ) == norient-1 )
            ix = ceil(posT.Mx(iF)*posScale / stepx);
            iy = ceil(posT.My(iF)*posScale / stepy);
            io = ceil( (posT.Mi(iF)+1) / stepo );
            found = true;
            weights(ix,iy,io) = weights(ix,iy,io) - negT.Mm(iF);
            break;
        end
    end
    if found
        continue;
    end
    ix = ceil( negT.Mx(iF)*negScale / stepx );
    iy = ceil( negT.My(iF)*negScale / stepy );
    io = ceil( (negT.Mi(iF)+1) / stepo );
    weights(ix,iy,io) = weights(ix,iy,io) - negT.Mm(iF);
    if detailed_xyo(ix,iy,io,1) == 0
        detailed_xyo(ix,iy,io,1:3) = [ceil(negT.Mx(iF)*negScale) ceil(negT.My(iF)*negScale) negT.Mi(iF)];
    end
end

posSym = zeros(120,120);
negSym = zeros(120,120);
commonSym = zeros(120,120);
for ix = 1:120/stepx
    for iy = 1:120/stepy
        for io = 1:norient/stepo
            
            if detailed_xyo(ix,iy,io,3) == 0
                continue;
            end
            val = weights(ix,iy,io);
            symbol = allsymbol{ detailed_xyo(ix,iy,io,3) };
            if val > weight_thres
                posSym(detailed_xyo(ix,iy,io,1)+(-h:h),detailed_xyo(ix,iy,io,2)+(-h:h)) = ...
                    max(posSym(detailed_xyo(ix,iy,io,1)+(-h:h),detailed_xyo(ix,iy,io,2)+(-h:h)) ...
                    , 1 * symbol);
            elseif val < -weight_thres
                negSym(detailed_xyo(ix,iy,io,1)+(-h:h),detailed_xyo(ix,iy,io,2)+(-h:h)) = ...
                    max(negSym(detailed_xyo(ix,iy,io,1)+(-h:h),detailed_xyo(ix,iy,io,2)+(-h:h)) ...
                    , 1 * symbol);
            else
                commonSym(detailed_xyo(ix,iy,io,1)+(-h:h),detailed_xyo(ix,iy,io,2)+(-h:h)) = ...
                    max(commonSym(detailed_xyo(ix,iy,io,1)+(-h:h),detailed_xyo(ix,iy,io,2)+(-h:h)) ...
                    , 1 * symbol);
            end
        end
    end
end

basisSym = zeros(120,120,3);
basisSym(:,:,2) = basisSym(:,:,2) - commonSym ;
%basisSym(:,:,2) = basisSym(:,:,2) - commonSym ;
%basisSym(:,:,3) = basisSym(:,:,3) - commonSym ;
basisSym(:,:,1) = basisSym(:,:,1) - posSym ;
% basisSym(:,:,3) = basisSym(:,:,3) - posSym ;
basisSym(:,:,3) = basisSym(:,:,3) - negSym ;
% basisSym(:,:,1) = basisSym(:,:,1) - 1.2*negSym ;
basisSym = 255 * (basisSym-min(basisSym(:)))/(max(basisSym(:))-min(basisSym(:)));
basisSym = uint8(basisSym);


%%
function basisSym = get_SVM_template(detailed_xyo,pos,neg,pos_thres,neg_thres,posSx,negSx,weight_thres,allfilter,allsymbol)
sx = 120; norient = 16;
examplefile = sprintf('svmtrain_%s_%s.txt',pos,neg);

	selectedx = [];
	selectedy = [];
	selectedo = [];

for j1 = 1:size(detailed_xyo,1)
	for j2 = 1:size(detailed_xyo,2)
		for j3 = 1:size(detailed_xyo,3)
			if detailed_xyo(j1,j2,j3,1) == 0
				continue;
			end
			xyo = detailed_xyo(j1,j2,j3,1:3);
		
			selectedx = [selectedx;xyo(1)];
			selectedy = [selectedy;xyo(2)];
			selectedo = [selectedo;xyo(3)];
		end
	end
end

if ~exist(examplefile,'file')

	fid = fopen( examplefile, 'w' );
	% compute M1 map
	folder = pos;
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
		if size(im,3) == 3
			im = rgb2gray(im);
		end
		I{i} = double(imresize(im,[sx,sx]));
	end
	allfiltered = applyfilter(I, allfilter);  % filter training images
	% -----------------------------------------------------------------
	% WHITENING
	% -----------------------------------------------------------------
	Upperbound = 6;
	Ctransform(nimage, norient, allfiltered, ones(1,nimage)*sx, ones(1,nimage)*sx, Upperbound);
	for j = 1:numel(allfiltered)
		allfiltered{j} = single(allfiltered{j});
	end
	subsampleM1 = 1;
	% prepare the svm input

	for i = 1:nimage
		count = 1;
		fprintf(fid, '%d ', 1);
		M1Map = mexc_ComputeMAX1( norient, allfiltered(i,:), 6/17, ...
		  	1, 17, subsampleM1 );
		for j1 = 1:size(detailed_xyo,1)
			for j2 = 1:size(detailed_xyo,2)
				for j3 = 1:size(detailed_xyo,3)
					if detailed_xyo(j1,j2,j3,1) == 0
						continue;
					end
					%disp( sprintf('%d %d %d', xyo(1), xyo(2), xyo(3)) );
					%drawnow;
					xyo = detailed_xyo(j1,j2,j3,1:3);
					val = M1Map{xyo(3)+1}( xyo(1), xyo(2) );
					fprintf(fid, '%d:%.5f ', count, val );
					count = count + 1;
				end
			end
		end
		fprintf(fid, '\r\n' );
	end

	% compute M1 map
	folder = neg;
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
		if size(im,3) == 3
			im = rgb2gray(im);
		end
		I{i} = double(imresize(im,[sx,sx]));
	end
	allfiltered = applyfilter(I, allfilter);  % filter training images
	% -----------------------------------------------------------------
	% WHITENING
	% -----------------------------------------------------------------
	Upperbound = 6;
	Ctransform(nimage, norient, allfiltered, ones(1,nimage)*sx, ones(1,nimage)*sx, Upperbound);
	for j = 1:numel(allfiltered)
		allfiltered{j} = single(allfiltered{j});
	end
	subsampleM1 = 1;
		  
	% prepare the svm input

	for i = 1:nimage
		count = 1;
		M1Map = mexc_ComputeMAX1( norient, allfiltered(i,:), 6/17, ...
		  	1, 17, subsampleM1 );
		fprintf(fid, '%d ', -1);
		for j1 = 1:size(detailed_xyo,1)
			for j2 = 1:size(detailed_xyo,2)
				for j3 = 1:size(detailed_xyo,3)
					if detailed_xyo(j1,j2,j3,1) == 0
						continue;
					end
					xyo = detailed_xyo(j1,j2,j3,1:3);
				
					val = M1Map{xyo(3)+1}( xyo(1), xyo(2) );
					fprintf(fid, '%d:%.5f ', count, val );
					count = count + 1;
				end
			end
		end
		fprintf(fid, '\r\n' );
	end

	fclose(fid);

end

modelfile = sprintf('svmmodel_%s_%s.txt', pos,neg);
if ~exist(modelfile,'file')
	system( sprintf('svm_learn -c 1 %s %s',examplefile,modelfile) );
	pause(0.5);
end


% load svm model
fid = fopen(modelfile, 'r');
		
% read line 10: number of SV
for line = 1:9
	tline = fgetl(fid);
end
nSV = fscanf(fid, '%d', 1) - 1;
disp(sprintf('nSV=%d',nSV));
tline = fgetl(fid);
tline = fgetl(fid);

D = length(selectedx);

lambdas = zeros(1,D);

% read line 12~end: SV's
for jSV = 1:nSV
	% read alpha * y
	alpha = fscanf( fid, '%f', 1 );
	vec = zeros(1,D);
	for d = 1:D
		feaID = fscanf( fid, '%d', 1 );
		colon = fscanf( fid, '%c', 1 );
		feaVal = fscanf( fid, '%f', 1 );
		vec(d) = feaVal;
	end
	lambdas = lambdas + vec * alpha;
	tline = fgetl(fid);
end


%
% visualize the svm template
%
indPos = find(lambdas>weight_thres);
symPos = displayGaborTemplate([120,120], selectedx(indPos), selectedy(indPos), selectedo(indPos), ...
	ones(1,length(indPos)), 0 );
indNeg = find(lambdas<-weight_thres);
symNeg = displayGaborTemplate([120,120], selectedx(indNeg), selectedy(indNeg), selectedo(indNeg), ...
	ones(1,length(indNeg)), 0 );
indCommon = find(lambdas>=-weight_thres & lambdas<=weight_thres);
symCommon = displayGaborTemplate([120,120], selectedx(indCommon), selectedy(indCommon), selectedo(indCommon), ...
	ones(1,length(indCommon)), 0 );
symCommon = double(symCommon);
symPos = double(symPos);
symNeg = double(symNeg);

basisSym = zeros(120,120,3);
basisSym(:,:,2) = basisSym(:,:,2) + symCommon;
basisSym(:,:,1) = basisSym(:,:,1) + symPos;
basisSym(:,:,3) = basisSym(:,:,3) + symNeg ;
basisSym = 255 * (basisSym-min(basisSym(:)))/(max(basisSym(:))-min(basisSym(:)));
basisSym = uint8(basisSym);

