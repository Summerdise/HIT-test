function [learnedModel] = learnOriHist_minVar(allfiltered,nF,h,backgroundIm,ax,ay)

% nF: maximum number of features to select

get_natural_statististics(ax);
load storedHI

% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_thres = 0.01;
nimage = size(allfiltered,1);
norient = size(allfiltered,2);
Upperbound = 6;
sx = size(allfiltered{1,1},1);
sy = size(allfiltered{1,1},2);

hI = cell(size(allfiltered,1),size(allfiltered,2));
for i = 1:size(hI,1)
    for j = 1:size(hI,2)
        hI{i,j} = zeros(size(allfiltered{1,1}));
    end
end

varMap = 1.5*ones(size(allfiltered{1,1}));

% compute map of orientation histogram
% mex CGraHist.c
sub = max(floor(ax/2),floor(h/2));
CGraHist(nimage,norient,allfiltered,hI,varMap,Upperbound,sx,sy,ax,ay,h,sub);
towrite = varMap;
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
imwrite(towrite,'varMap.jpg');
varMapSaved = varMap;

% pursuit the least-variance positions
inhibited = zeros(sx,sy);
selectedMeanHist = zeros(nF,norient);
maxValue = max(varMap(:));
minVar = zeros(nF,1);
Mx = zeros(nF,1); My = zeros(nF,1);

for k = 1:nF
    [minVar(k) ind] = min(varMap(:));
    
    col = ceil(ind/sx); My(k) = col;
    row = ind - (col-1)*sx; Mx(k) = row;
    % calculate the mean gradient histogram
    hists = zeros(nimage,norient);
    for o = 1:norient
        for j = 1:nimage
            hists(j,o) = hI{j,o}(row,col);
        end
    end
    selectedMeanHist(k,:) = sum(hists,1)/nimage;
    
    % re-estimate the normalizing constant
    responses = sum((storedHI - repmat(selectedMeanHist(k,:),[size(storedHI,1),1])).^2,2);
    dZ(k) = mean(exp( - responses / 2 / minVar(k) ) / (sqrt(2*pi*minVar(k)))); % additional normalizing constant
    
    if -1/2-log( sqrt(2*pi.*minVar(k)) )-log(dZ(k)) < gain_thres || minVar(k) > maxValue
        nF = k-1;
        break;
    end
    
    % do inhibition
    left = max(1,col-ay);
    right = min(sy,col+ay);
    top = max(1,row-ax);
    bottom = min(sx,row+ax);
    inhibited(top:bottom,left:right) = 1;
    varMap(top:bottom,left:right) = maxValue + 1;
end
minVar = minVar(1:nF); Mx = Mx(1:nF); My = My(1:nF); selectedMeanHist = selectedMeanHist(1:nF,:); 
gain = (-1/2-log(sqrt(2*pi.*minVar))) - log(dZ(1:nF)');

towrite = double(inhibited) .* double(backgroundIm);
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));

learnedModel = struct('var',minVar,'meanHist',selectedMeanHist,...
    'inhibited',inhibited,'varMap',varMapSaved,'sym',towrite,'Mx',Mx,'My',My...
    , 'ax',ax,'ay',ay,'sx',sx,'sy',sy,'gain',gain,'dZ',dZ);

