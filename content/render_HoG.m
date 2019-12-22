% render HoG feature maps
clear
close all
addpath phog
load config1

bin = 8;
angle = 180;
L = 4;
% starting = 169; % for L = 3
starting = 681;
% tSize = [360 360];
tSize = 1;

folder = 'wolfHead';
Iname = dir(sprintf('%s/*.jpg',folder));
collected = [];
for imageID = 1:length(Iname)
	im = imread( sprintf( '%s/%s', folder, Iname(imageID).name ) );
    im = imresize( im, tSize );
	roi = [1;size(im,1);1;size(im,2)];
	HoG = anna_phog(sprintf( '%s/%s', folder, Iname(imageID).name ),bin,angle,L,roi);
	partialHoG = HoG(starting:end);
    if numel(partialHoG)/bin ~= (2^L)^2
        continue;
    end
	collected = [collected;partialHoG(:)'];
	
end
partialHoG = mean( collected, 1 );
partialHoG = reshape(partialHoG,bin,numel(partialHoG)/bin);
sideWidth = sqrt( size(partialHoG,2) );
canvas = zeros( sideWidth*(2*h+1) );
for row0 = 1:sideWidth
	for col0 = 1:sideWidth
		rows = (row0-1)*(2*h+1) + ( 1:(2*h+1) );		
		cols = (col0-1)*(2*h+1) + ( 1:(2*h+1) );
		vec = partialHoG(:,row0+(col0-1)*sideWidth);
		vec = vec / (sum(vec)+1e-10);
		for o = 1:bin
			o2 = bin - o + 1;
			canvas( rows, cols ) = canvas( rows, cols ) + vec(o) * allsymbol{o2*2-1} * sum(partialHoG(:,row0+(col0-1)*sideWidth));
		end
	end
end
towrite = -canvas;
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
imwrite( towrite, sprintf('hog_%s.png',folder) );



folder = 'hedgehog';
Iname = dir(sprintf('%s/*.jpg',folder));
collected = [];
for imageID = 1:length(Iname)
	im = imread( sprintf( '%s/%s', folder, Iname(imageID).name ) );
    im = imresize( im, tSize );
	roi = [1;size(im,1);1;size(im,2)];
	HoG = anna_phog(sprintf( '%s/%s', folder, Iname(imageID).name ),bin,angle,L,roi);
	partialHoG = HoG(starting:end);
    if numel(partialHoG)/bin ~= (2^L)^2
        continue;
    end
	collected = [collected;partialHoG(:)'];
	
end
partialHoG = mean( collected, 1 );
partialHoG = reshape(partialHoG,bin,numel(partialHoG)/bin);
sideWidth = sqrt( size(partialHoG,2) );
canvas = zeros( sideWidth*(2*h+1) );
for row0 = 1:sideWidth
	for col0 = 1:sideWidth
		rows = (row0-1)*(2*h+1) + ( 1:(2*h+1) );
		cols = (col0-1)*(2*h+1) + ( 1:(2*h+1) );
		vec = partialHoG(:,row0+(col0-1)*sideWidth);
		vec = vec / (sum(vec)+1e-10);
		for o = 1:bin
			o2 = bin - o + 1;
			canvas( rows, cols ) = canvas( rows, cols ) + vec(o) * allsymbol{o2*2-1} * sum(partialHoG(:,row0+(col0-1)*sideWidth));
		end
	end
end
towrite = -canvas;
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
imwrite( towrite, sprintf('hog_%s.png',folder) );



folder = 'pigeonHead';
Iname = dir(sprintf('%s/*.jpg',folder));
collected = [];
for imageID = 1:length(Iname)
	im = imread( sprintf( '%s/%s', folder, Iname(imageID).name ) );
    im = imresize( im, tSize );
	roi = [1;size(im,1);1;size(im,2)];
	HoG = anna_phog(sprintf( '%s/%s', folder, Iname(imageID).name ),bin,angle,L,roi);
	partialHoG = HoG(starting:end);
    if numel(partialHoG)/bin ~= (2^L)^2
        continue;
    end
	collected = [collected;partialHoG(:)'];
	
end
partialHoG = mean( collected, 1 );
partialHoG = reshape(partialHoG,bin,numel(partialHoG)/bin);
sideWidth = sqrt( size(partialHoG,2) );
canvas = zeros( sideWidth*(2*h+1) );
for row0 = 1:sideWidth
	for col0 = 1:sideWidth
		rows = (row0-1)*(2*h+1) + ( 1:(2*h+1) );
		cols = (col0-1)*(2*h+1) + ( 1:(2*h+1) );
		vec = partialHoG(:,row0+(col0-1)*sideWidth);
		vec = vec / (sum(vec)+1e-10);
		for o = 1:bin
			o2 = bin - o + 1;
			canvas( rows, cols ) = canvas( rows, cols ) + vec(o) * allsymbol{o2*2-1} * sum(partialHoG(:,row0+(col0-1)*sideWidth));
		end
	end
end
towrite = -canvas;
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
imwrite( towrite, sprintf('hog_%s.png',folder) );



folder = 'pigHead';
Iname = dir(sprintf('%s/*.jpg',folder));
collected = [];
for imageID = 1:length(Iname)
	im = imread( sprintf( '%s/%s', folder, Iname(imageID).name ) );
    im = imresize( im, tSize );
	roi = [1;size(im,1);1;size(im,2)];
	HoG = anna_phog(sprintf( '%s/%s', folder, Iname(imageID).name ),bin,angle,L,roi);
	partialHoG = HoG(starting:end);
    if numel(partialHoG)/bin ~= (2^L)^2
        continue;
    end
	collected = [collected;partialHoG(:)'];
	
end
partialHoG = mean( collected, 1 );
partialHoG = reshape(partialHoG,bin,numel(partialHoG)/bin);
sideWidth = sqrt( size(partialHoG,2) );
canvas = zeros( sideWidth*(2*h+1) );
for row0 = 1:sideWidth
	for col0 = 1:sideWidth
		rows = (row0-1)*(2*h+1) + ( 1:(2*h+1) );
		cols = (col0-1)*(2*h+1) + ( 1:(2*h+1) );
		vec = partialHoG(:,row0+(col0-1)*sideWidth);
		vec = vec / (sum(vec)+1e-10);
		for o = 1:bin
			o2 = bin - o + 1;
			canvas( rows, cols ) = canvas( rows, cols ) + vec(o) * allsymbol{o2*2-1} * sum(partialHoG(:,row0+(col0-1)*sideWidth));
		end
	end
end
towrite = -canvas;
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
imwrite( towrite, sprintf('hog_%s.png',folder) );
