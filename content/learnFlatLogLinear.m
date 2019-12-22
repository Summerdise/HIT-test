function [learnedModel] = learnFlatLogLinear(LocalFlatLearn,maxNumFeature,ax)

load nlf_flat
storedMean = e;
storedLambda = lam;
storedLogZ = lz;

% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_thres = 0.5;
nimage = size(LocalFlatLearn,1);
sx = size(LocalFlatLearn{1,1},1);
sy = size(LocalFlatLearn{1,1},2);

% compute map of orientation histogram

% compute the map of information gain
inhibited = zeros(sx,sy);
templateSize = [sx sy];

meanLocalFlatMap = zeros(sx, sy);
for i = 1:nimage
	meanLocalFlatMap(:,:) = meanLocalFlatMap(:,:) + 1/nimage * LocalFlatLearn{i};
end

Mx = zeros(1, maxNumFeature);  % position of selected Gabor
My = zeros(1, maxNumFeature); 
Mm = zeros(1, maxNumFeature);  % mean
lambda = zeros(1, maxNumFeature);  % lambda
logZ = zeros(1, maxNumFeature);  % logZ
gain = zeros(maxNumFeature, 1);  % coding gain for selected Gabor

gainMap = zeros(templateSize(1),templateSize(2));
ind = zeros(templateSize(1),templateSize(2));
sub = 4;
for x = sub/2+1:sub:templateSize(1)-sub/2
	for y = sub/2+1:sub:templateSize(2)-sub/2
        [dif ind(x,y)] = min( abs(storedMean - meanLocalFlatMap(x,y)) );
        ind(x,y) = min(length(storedMean),max(1,ind(x,y)));
        gainMap(x,y) = storedMean(ind(x,y)) * storedLambda(ind(x,y)) - storedLogZ(ind(x,y));
	end
end


gainMapSaved = gainMap;

gainMap([1:ax*2 templateSize(1)-2*ax+1:end],:) = 0;
gainMap(:,[1:ax*2 templateSize(1)-2*ax+1:end]) = 0;

% pursuit the largest-gain positions

nF = maxNumFeature;
for k = 1:nF
	[maxGain here] = max(gainMap(:) .* double(1-inhibited(:))  );
	binNo = ind(here);

	if maxGain <= gain_thres
		nF = k-1;
		break
	end

	col = ceil(here/templateSize(1));
	My(k) = col;
	row = here - (col-1)*templateSize(1); 
	Mx(k) = row;
	
	gain(k) = maxGain;
	Mm(k) = storedMean(binNo);
	lambda(k) = storedLambda(binNo);
	logZ(k) = storedLogZ(binNo);
	
	% do inhibition
	left = max(1,col-ax);
	right = min(templateSize(2),col+ax);
	top = max(1,row-ax);
	bottom = min(templateSize(1),row+ax);
    inhibited(top:bottom,left:right) = 1;
    gainMap(top:bottom,left:right) = 0;
    disp(sprintf('%d %d %.3f %.3f',Mx(k),My(k),Mm(k),gain(k)));
end
Mm = Mm(1:nF); Mx = Mx(1:nF); My = My(1:nF);
gain = gain(1:nF);
lambda = lambda(1:nF);
logZ = logZ(1:nF);

% save results

learnedModel = struct(...
    'inhibited',inhibited,'gainMap',gainMapSaved,'Mx',Mx,'My',My...
    ,'ax',ax,'sx',sx,'sy',sy,'gain',gain,'logZ',logZ,'lambda',lambda,'Mm',Mm,...
    'meanLocalFlatMap',meanLocalFlatMap);

