function [allfiltered allfilteredr allfilteredi] = applyfilter(I, allfilter)
% filter image by a bank of filters
% I: input images
% allfilter: filter bank


nimage = size(I, 2);    % number of training images
[sx, sy] = size(I{1});  % size of images
norient = size(allfilter, 2);     % number of orientations
h = (size(allfilter{1}, 1)-1)/2;  % half size of filters
allfiltered = cell(nimage, norient);  % filtered images
allfilteredr = cell(nimage, norient);
allfilteredi = cell(nimage, norient);
for i = 1 : nimage
   tot = 0.; 
   for o = 1 : norient
       filtered = filter2(allfilter{1, o}, I{i}, 'same');
       re = real(filtered); im = imag(filtered);
       allfilteredr{i,o} = re; allfilteredi{i,o} = im;
       energy = re.*re + im.*im;  % compute the local energy
%        energy([1:h sx-h:sx],:) = 0; energy(:,[1:h sy-h:sy]) = 0; % modified !!!
       allfiltered{i, o} = energy;
       tot = tot + sum(sum(energy((h+1):(sx-h-1), (h+1):(sy-h-1))))/(sx-2*h-1)/(sy-2*h-1); 
   end
   ave = tot/norient; 
   for o = 1 : norient
       allfiltered{i, o} = allfiltered{i, o}/ave;  % normalizing by whitening
   end
end

