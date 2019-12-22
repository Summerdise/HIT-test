function [allfiltered] = applyfilterfftsame(I, allfilter)
% filter image by a bank of filters
% I: input images
% allfilter: filter bank
% threshold: filter response is set to zero if below threshold

nimage = size(I, 2);    % number of training images
norient = size(allfilter, 2);     % number of orientations
h = (size(allfilter{1}, 1)-1)/2;  % half size of filters
allfiltered = cell(nimage, norient);  % filtered images
for (i = 1 : nimage)
    [sx sy] = size(I{i}); % size of images
   fftI{i} = fft2(I{i}, sx+h+h, sy+h+h);
end
for i = 1:nimage
    [sx sy] = size(I{i}); % size of images
    for  (o = 1 : norient)
       fftf{i,o} = fft2(allfilter{o}, sx+h+h, sy+h+h); 
    end
end

for (i = 1 : nimage)
    [sx sy] = size(I{i}); % size of images
   tot = 0.;
   fftIi = fftI{i}; 
   for (o = 1 : norient)
      fftfo = fftf{i,o}; 
      out = ifft2(fftIi.*fftfo);
      filtered = out(h+1:h+sx, h+1:h+sy); 
       re = real(filtered); im = imag(filtered); 
       energy = re.*re + im.*im;  % compute the local energy
       energy([1:h sx-h:sx],:) = 0; energy(:,[1:h sy-h:sy]) = 0; % modified !!!
       allfiltered{i, o} = energy; 
       tot = tot + sum(sum(energy((h+1):(sx-h-1), (h+1):(sy-h-1))))/(sx-2*h-1)/(sy-2*h-1); 
   end
   ave = tot/norient; 
   for o = 1 : norient
       allfiltered{i, o} = allfiltered{i, o}/ave;  % normalizing by whitening
   end
end

% -----------------------------------------------------------------
% LOCAL ADAPTIVE WHITENING
% -----------------------------------------------------------------
% mex CLocalAdapt.c
% allfiltered2 = cell(size(allfiltered,1),size(allfiltered,2));
% for i = 1:size(allfiltered,1)
%     for j = 1:size(allfiltered,2)
%         allfiltered2{i,j} = allfiltered{i,j} + 0;
%     end
% end
% CLocalAdapt(nimage,norient,allfiltered2,allfiltered,sx,sy);

% for i = 1 : nimage
%    tot = zeros(sx,sy); 
%    fftIi = fftI{i};
%    for o = 1 : norient
%        fftfo = fftf{o};
%        out = ifft2(fftIi.*fftfo);
%        filtered = out(h+1:h+sx, h+1:h+sy);
%        % filtered = filter2(allfilter{1, o}, I{i}, 'same');
%        re = real(filtered); im = imag(filtered); 
%        energy = re.*re + im.*im;  % compute the local energy
%        allfiltered{i, o} = energy; 
%        tot = tot + energy(:, :); 
%    end
%    ave = tot/norient + 1e-6;
%    for o = 1 : norient
%        allfiltered{i, o} = allfiltered{i, o}./ave;  % normalizing by whitening
%    end
% end

