% configuration file for matlab
clear
scale = 0.7;
norient = 16;  % number of orientations

T = 0.;        % threshold for responses 
sub = 2;    % subsampling rate for Gabor
Lrange = 3; % allowed shifting in location (times sub pixels)
Orange = 1; % allowed shifting in orientation
%%% Note: Lrange*sub is the pixels allowed to shift along normal
%%%       Lrange*sub should be the same in learning and testing!
[allfilter allsymbol] = makefilter(scale, norient);  % generate Gabor filters
h = (size(allfilter{1}, 1)-1)/2;  % half size of Gabor
Upperbound = 6;
save('config1.mat');
