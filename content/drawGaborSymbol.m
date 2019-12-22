function im = drawGaborSymbol(im, allsymbol, row, col, orientationIndex, intensity)

nRow = size(im,1);
nCol = size(im,2);
h = (size(allsymbol{1}, 1)-1)/2;  % half size of Gabor

% position in image
top = row-h;
down = row+h;
left = col-h;
right = col+h;

% top2, down2, ... is the position in filter

% deal with out-of-boundary cases
if top < 1
    top2 = 1-top + 1;
    top = 1;
else
    top2 = 1;
end

if down > nRow
    down2 = down-nRow;
    down = nRow;
else
    down2 = size(allsymbol{1},1);
end

if left < 1
    left2 = 1-left + 1;
    left = 1;
else
    left2 = 1;
end

if right > nCol
    right2 = right - nCol;
    right = nCol;
else
    right2 = size(allsymbol{1},1);
end

to_draw = allsymbol{orientationIndex} * intensity;
mask = im(top:down,left:right) >= to_draw(top2:down2,left2:right2);
im(top:down,left:right) = (1-mask).*to_draw(top2:down2,left2:right2) + mask.*im(top:down,left:right);

