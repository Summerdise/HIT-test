function towrite = displayGaborTemplate(templateSize,Mx,My,Mi,Mm,normalize)
if nargin < 6
	normalize = true;
end

sx = templateSize(1); sy = templateSize(2);

scale = .7;  % scale of Gabors
norient = 16;  % number of orientations
[allfilter, allsymbol] = makefilter(scale, norient);  % generate Gabor filters 

sym = zeros(sx,sy);
nF = length(Mx);
for k = 1:nF
    ori = Mi(k) + 1;
    col = My(k);
    row = Mx(k);
    sym = drawGaborSymbol( sym, allsymbol, row, col, ori, sqrt(Mm(k)) );
end

towrite = -sym;
if normalize
	towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
end
