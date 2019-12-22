function learnedModel = ABlearn2(dataWeight,allfiltered,allsymbol,allfilter,T,Lrange,Orange,sub,Totalsketch)

load('nlf.mat'); 
epsilon = .1;  % allowed correlation between selected Gabors
h = (size(allfilter{1}, 1)-1)/2;  % half size of Gabor
C = corr(allfilter, epsilon);    % generate the inhibition maps

%%% Note: Lrange*sub is the pixels allowed to shift along normal
%%%       Lrange*sub should be the same in learning and testing!
SHUTUP = 0.;
% SHUTUP is to neglect the positions with low average responses
% VARIABLES FOR MEX-C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mi = zeros(1, Totalsketch);  % orientation of selected Gabor
Mx = zeros(1, Totalsketch);  % position of selected Gabor
My = zeros(1, Totalsketch); 
Mm = zeros(1, Totalsketch);  % lambda
Mm1 = zeros(1, Totalsketch);  % logZ
sx = size(allfiltered{1,1},1);
sy = size(allfiltered{1,1},2);
sym = zeros(sx,sy);  % symbolic plot of selected Gabor
nimage = size(allfiltered,1);
norient = size(allfiltered,2);
Asym = cell(1, nimage);   % symbolic plot for each image with active Gabor
for i = 1 : nimage
    Asym{i} = zeros(sx,sy);  
end
tss = zeros(nimage, 1);  % total scores for testing images
gain = zeros(Totalsketch, 1);  % coding gain for selected Gabor

binnum = 50; % number of bins 
% MEX-C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mex Clearn2.c;   % compile C code
Upperbound = 6;
Clearn2(nimage, norient, dataWeight, allfiltered, C, h, sx, sy, T, Mi, Mx, My, Mm, Mm1, allsymbol(1, :), sym, Asym, Lrange, Orange, sub, Totalsketch, Upperbound, tss, gain, SHUTUP, binnum, M, lam, e, lz);

% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
learnedModel = struct('Mm',Mm,'Mm1',Mm1,'Mi',Mi,'Mx',Mx,'My',My,...
    'Lrange',Lrange,'Orange',Orange,'h',h,'sub',sub,'sx',sx,'sy',sy,...
    'sym',sym,'T',T,'Upperbound',Upperbound,'epsilon',epsilon,...
    'Totalsketch',Totalsketch,'gain',gain);
