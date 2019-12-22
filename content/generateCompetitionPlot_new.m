function generateCompetitionPlot_new(folder,sxBySy,sxBySyTex)
if nargin < 3
    sxBySyTex = 150 * 150;
end

gain_thres = 0.1;
nChannel = 0;

load(sprintf('rawmodel_basis_%s_size%d.mat',folder,sxBySy),'template'); nChannel = nChannel + 1; 
template_all{nChannel} = template;

ax = 10;
load(sprintf('rawmodel_tex_%s_size%d_ax%d.mat',folder,sxBySyTex,ax),'template'); nChannel = nChannel + 1; 
template_all{nChannel} = template;

ax = 10;
load(sprintf('rawmodel_flat_%s_size%d_ax%d.mat',folder,sxBySyTex,ax),'template'); nChannel = nChannel + 1; 
template_all{nChannel} = template;

ax = 20;
isColor = true;
filename = sprintf('rawmodel_color_%s_size%d_ax%d.mat',folder,sxBySyTex,ax);
if ~exist(filename,'file')
	isColor = false;
end
if isColor
	load(filename,'template'); nChannel = nChannel + 1; 
	template_all{nChannel} = template;
end

gainBasis = template_all{1}.gain; gainBasis = gainBasis(gainBasis>gain_thres);
gainTex = template_all{2}.gain; gainTex = gainTex (gainTex >gain_thres);
gainFlat = template_all{3}.gain; gainFlat = gainFlat (gainFlat >gain_thres);
gainvec = [gainBasis; gainTex; gainFlat];
if isColor
	gainColor = template_all{4}.gain; gainColor = gainColor (gainColor >gain_thres);
	gainvec = [gainvec; gainColor ];
end
typePerFeature = [ones(1,length(gainBasis)), 2 * ones(1,length(gainTex )) , 3 * ones(1,length(gainFlat))];
if isColor
	typePerFeature = [typePerFeature, 4 * ones(1,length(gainColor)) ]';
end

gainMat = [typePerFeature(:) gainvec];
[sorted ind] = sortrows(gainMat,-2);
h=figure;hold on
if sum(sorted(:,1)==1)>0
    bar(1:size(sorted,1),sorted(:,2).*(sorted(:,1)==1),'w');
end
if sum(sorted(:,1)==2)>0
    bar(1:length(ind),sorted(:,2).*(sorted(:,1)==2),'r'); 
end
if sum(sorted(:,1)==3)>0
    bar(1:length(ind),sorted(:,2).*(sorted(:,1)==3),'b'); 
end
if isColor && sum(sorted(:,1)==4)>0
    bar(1:length(ind),sorted(:,2).*(sorted(:,1)==4),'g'); 
end

xlim([0 40]);
set(gcf, 'PaperPositionMode', 'auto');
if isColor
	legend({'Sketch','Texture','Flat','Color'},'FontSize',30);
else
	legend({'Sketch','Texture','Flat'},'FontSize',30);
end
saveas(h,['competition/' folder '.png']);
disp('press any key to continue'); pause(1);
close(h);
