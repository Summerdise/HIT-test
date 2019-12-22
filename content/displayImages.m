function im = displayImages(images,ncol,bx,by,normalize)
% displayImage - display a set of images of the same size (bx, by) on 
%   a specified tabular.
%
% ncol: number of columns

if nargin < 5
    normalize = true;
end

n = length(images);
nrow = ceil(n/ncol);
widthMargin = 5;
width = by*ncol + (ncol-1)*widthMargin;
heightMargin = 5;
height = bx*nrow + (nrow-1)*heightMargin;

im = 255-zeros(height,width,3);
for i=1:n
    row = ceil(i/ncol);
    col = i-(row-1)*ncol;
    startx = (row-1)*(bx+heightMargin);
    starty = (col-1)*(by+widthMargin);
    if isempty(images{i})
    	continue;
    end
    towrite = imresize(double(images{i}),[bx,by]);
    if normalize
        towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
    else
        towrite = uint8(towrite);
    end
    if size(towrite,3) == 1
    	im((startx+1):(startx+bx),(starty+1):(starty+by),:) = repmat(towrite,[1,1,3]);
    else
    	im((startx+1):(startx+bx),(starty+1):(starty+by),:) = towrite;
    end
end
im = uint8(im);
