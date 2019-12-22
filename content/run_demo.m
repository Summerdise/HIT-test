clear
close all


mex Clearn2.c
mex Ctransform.c
mex Chistog.c
mex CGraHist.c

nrow = 3; ncol = 2;
%%

handle = figure; hold on
text(20,20,'Please press Enter to continue');
count = 1;

folder = 'bearHead';
images = get_examples(folder);
for i = 1:4
    subplot(nrow,ncol,count);
    imshow(images{i},[]);
    count = count + 1;
end
subplot(nrow,ncol,count);
disp(['learning mixed template for ' folder]);
sxBySy = 70*70; sxBySyTex = 120*120;
sk_thres = 0.6;
tex_thres = 1;
flat_thres = 0.8;
color_thres = 0.6;
learn_mixed(folder,sxBySy,sk_thres,tex_thres,flat_thres,color_thres,sxBySyTex);
generateCompetitionPlot_new(folder,sxBySy,sxBySyTex);

clear

%}
