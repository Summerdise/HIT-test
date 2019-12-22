function learnedModel = learnOriHist_table(allfiltered,nF,h,backgroundIm,ax,ay)

% nF: number of features to select
% we do not use dataWeight for now (to modify!!)
 
% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 
varMap = 0.05*ones(size(allfiltered{1,1}));
 
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
    if minVar(k) > 0.05 - 1e-3 % -1/2-log(sqrt(2*pi.*minVar(k))) < 0.1 ||
        nF = k-1;
        break;
    end
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
    % do inhibition
    left = max(1,col-ay);
    right = min(sy,col+ay);
    top = max(1,row-ax);
    bottom = min(sx,row+ax);
    inhibited(top:bottom,left:right) = 1;
    varMap(top:bottom,left:right) = maxValue + 1;
end
minVar = minVar(1:nF); Mx = Mx(1:nF); My = My(1:nF); selectedMeanHist = selectedMeanHist(1:nF,:); 

% prepare |h-h*|
hMetric = cell(1,nF); qMetric = cell(1,nF);
load(['oriHistDB_ax' num2str(ax) '.mat'], 'oriHistDB');
bins_metric = 0:0.001:1;
for iF = 1:nF
    hMetricVec = [];
    for i = 1:min(size(oriHistDB,2),50)
        met = sum( (oriHistDB{i} - repmat(selectedMeanHist(iF,:),[size(oriHistDB{i},1) 1])).^2,2);
        hMetricVec = [hMetricVec;met];
    end
    
    hMetric{iF} = hMetricVec;
    qMetric{iF} = hist(hMetricVec,bins_metric);
    qMetric{iF} = qMetric{iF}/sum(qMetric{iF});
end
% 
lambda = zeros(nF,1); logZ = zeros(nF,1); gain = zeros(nF,1);
lambda_range = -100 : 1 : -1; Z = zeros(1,length(lambda_range));
e = zeros(1,length(lambda_range));
for iF = 1:nF
    for j = 1:length(lambda_range)
        cand_lambda = lambda_range(j);
        p = exp(cand_lambda*bins_metric).*qMetric{iF}; 
        Z(j) = sum(p); p = p/Z(j);
%         Z(j) = sum(p*(bins_metric(2)-bins_metric(1))); p = p/Z(j);
        e(j) = sum(bins_metric.*p);
    end
    ind = sum(e<minVar(iF));
    if ind == length(e)
        lambda(iF) = lambda_range(end);
        logZ(iF) = log(Z(end));
         
    elseif ind == 0
        lambda(iF) = lambda_range(1);
        logZ(iF) = log(Z(1));

    else
        lambda(iF) = lambda_range(ind) + (lambda_range(ind+1) - lambda_range(ind)) * (minVar(iF)-e(ind));
        logZ(iF) = log(Z(ind)) + (log(Z(ind+1)) - log(Z(ind))) * (minVar(iF)-e(ind));

    end
    gain(iF) = lambda(iF) * minVar(iF) - logZ(iF);
end

[gain ind] = sort(gain,'descend');
Mx = Mx(ind); My = My(ind); logZ = logZ(ind); lambda = lambda(ind);
minVar = minVar(ind);
tmp = selectedMeanHist;
for iF = 1:nF
    selectedMeanHist(iF,:) = tmp(ind(iF),:);
end
 
towrite = double(inhibited) .* double(backgroundIm);
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
 
learnedModel = struct('var',minVar,'meanHist',selectedMeanHist,...
    'inhibited',inhibited,'varMap',varMapSaved,'sym',towrite,'Mx',Mx,'My',My...
    , 'ax',ax,'ay',ay,'sx',sx,'sy',sy,'gain',gain,'lambda',lambda,'logZ',logZ);

