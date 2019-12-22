function storedHI = get_natural_statististics(ax)

% store the histogram features extracted from random natural images
if exist('storedHI.mat','file')
   load storedHI;
   return;
end

load config1

negFolder = 'rand';  % image input folder
Iname = dir([negFolder '/*.jpg']);
nimage = min(length(Iname),200); % #images
Iname = Iname(1:nimage);
Isize = zeros(nimage, 2);
I = cell(1,nimage);
colorI = cell(1,nimage);
for i = 1 : nimage
    tmpIm = imread([negFolder '/' Iname(i).name]);
    
    if size(tmpIm,3) == 3
    	colorI{i} = tmpIm;
        tmpIm = rgb2gray(tmpIm);
    end
    I{i} = double(tmpIm);
    Isize(i, :) = size(I{i});
end
sx = min(Isize(:, 1)); sy = min(Isize(:, 2)); 
for i = 1 : nimage
    I{i} = I{i}(1:sx, 1:sy); 
end   % make the sizes of the testing images the same


% FILTERING IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allfilter = makefilter(scale, norient);  % generate Gabor filters 
h = (size(allfilter{1}, 1)-1)/2;  % half size of Gabor
allfiltered = applyfilterfftsame(I, allfilter);  % filter training images

Ctransform(nimage, norient, allfiltered, Isize(:,1), Isize(:,2), Upperbound);

hI = cell(size(allfiltered,1),size(allfiltered,2));
for i = 1:size(hI,1)
    for j = 1:size(hI,2)
        hI{i,j} = zeros(size(allfiltered{1,1}));
    end
end

varMap = 1.5*ones(size(allfiltered{1,1}));

% compute map of orientation histogram
sub = 30;
CGraHist(nimage,norient,allfiltered,hI,varMap,Upperbound,sx,sy,ax,ax,h,sub);

count = 1;
for i = 1:nimage
    for xx = 1+h+ax:sub:sx-ax-h
        for yy = 1+h+ax:sub:sy-ax-h
            tmp = zeros(norient,1);
            for o = 1:norient
                tmp(o) = hI{i,o}(xx,yy);
            end
            storedHI{count} = tmp;
            count = count + 1;
        end
    end
end

storedHI = cell2mat(storedHI);
storedHI = storedHI';

save('storedHI.mat','storedHI');


% compute map of local average RGB
localAveFilter = ones(20,20,'single')/single(20)/single(20);
scalar = 1.0 / 250;
count = 1;
storedHI = {};
for i = 1:nimage
	if isempty( colorI{i} )
		continue;
	end
	
	R = colorI{i}(:,:,1);
	G = colorI{i}(:,:,2);
	B = colorI{i}(:,:,3);
	R = filter2(localAveFilter,R,'same');
	G = filter2(localAveFilter,G,'same');
	B = filter2(localAveFilter,B,'same');

    for xx = 1+h+ax:sub:sx-ax-h
        for yy = 1+h+ax:sub:sy-ax-h
            storedHI{count} = [R(xx,yy); G(xx,yy); B(xx,yy)] * scalar;
            count = count + 1;
        end
    end
end
storedHI = cell2mat(storedHI);
storedHI = storedHI';
save('storedColorHI.mat','storedHI');
