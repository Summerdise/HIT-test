function learn_mixed(folder,sxBySy,sketch_thres,tex_thres,flat_thres,color_thres,sxBySyTex)

% ==================================================
% parameters
% ==================================================

if nargin < 7
    sxBySyTex = 150 * 150;
end


load config1 norient h allsymbol

txByTy = 250*250;

% ==================================================
% learn and illustrate the sketch template
% ==================================================

modelfile = sprintf('rawmodel_basis_%s_size%d.mat',folder,sxBySy);
if ~exist(modelfile,'file')
    tic
    [template bx by] = learnBasis(folder,sxBySy);
    disp(sprintf('selecting sketch features : %.3f seconds',toc));
else
    load(sprintf('rawmodel_basis_%s_size%d.mat',folder,sxBySy),'template','bx','by');
end

tx = floor(sqrt(txByTy/bx/by)*bx); ty = floor(txByTy/tx);

basisSym = zeros(bx,by);

Mx = template.Mx; My = template.My; Mo = template.Mi + 1; gain = template.gain;

for iF = 1:sum(gain>sketch_thres)
    val = template.Mm(iF);
    symbol = allsymbol{Mo(iF)};
    basisSym(Mx(iF)+(-h:h),My(iF)+(-h:h)) = ...
        max(basisSym(Mx(iF)+(-h:h),My(iF)+(-h:h)) ...
        , 1 * symbol);
end

% ==========================================================
% learn and illustrate the texture template
% ==========================================================
ax = 10; ay = ax;
modelfile = sprintf('rawmodel_tex_%s_size%d_ax%d.mat',folder,sxBySyTex,ax);
if ~exist(modelfile,'file')
    tic
    [template bx by] = learnTex(folder,sxBySyTex,ax);
    disp(sprintf('selecting texture features : %.3f seconds',toc));
else
    load(modelfile,'template','bx','by');
end
tx = bx; ty = by;
Mx = template.Mx; My = template.My; gain = template.gain;

texSym = zeros(bx,by);
for iF = 1:length(gain)
    if sum(gain(iF)<tex_thres)
        continue;
    end
    vec = template.meanHist(iF,:);
    vec = vec.^4 + 1e-10;
    vec = vec/sum(vec);
    
    small_symbol = zeros( h*2+1, h*2+1 );
    for ori = 1:norient
    	small_symbol = small_symbol + vec(ori) * allsymbol{ori};
    end
    
    paved_symbol = zeros( size(small_symbol)*2 + 1 );
    paved_symbol(1:2*h+1,1:2*h+1) = small_symbol;
    paved_symbol(2*h+1+(1:2*h+1),1:2*h+1) = small_symbol;
    paved_symbol(1:2*h+1,2*h+1+(1:2*h+1)) = small_symbol;
    paved_symbol(2*h+1+(1:2*h+1),2*h+1+(1:2*h+1)) = small_symbol;
    paved_symbol = imresize(paved_symbol,[ax+2*h+1 ax+2*h+1]);
    row_ind = ( Mx(iF) - floor(size(paved_symbol,1)/2) -1 ) + (1:size(paved_symbol,1));
    col_ind = ( My(iF) - floor(size(paved_symbol,2)/2) -1 ) + (1:size(paved_symbol,2));
    texSym(row_ind,col_ind) = ...
        max(texSym(row_ind,col_ind) ...
        , paved_symbol );

    
end


% ==========================================================
% learn and illustrate the flatness template
% ==========================================================
ax = 10; ay = ax;
sxBySyFlat = sxBySyTex;
modelfile = sprintf('rawmodel_flat_%s_size%d_ax%d.mat',folder,sxBySyTex,ax);
if ~exist(modelfile,'file')
    tic
    [template bx by] = learnFlat(folder,sxBySyFlat,ax);
    disp(sprintf('selecting flatness features : %.3f seconds',toc));
else
    load(modelfile,'template','bx','by');
end
tx = bx; ty = by;
Mx = template.Mx; My = template.My; gain = template.gain;

flatSym = zeros(bx,by);
for iF = 1:length(gain)
    if sum(gain(iF)<flat_thres)
        continue;
    end
    flatSym(Mx(iF)+(-ax:ax),My(iF)+(-ax:ax)) = 1;
end


% ==========================================================
% learn and illustrate the color template
% ==========================================================

ax = 20; ay = ax;
ax_color = ax;
sxBySyColor = sxBySyTex;
modelfile = sprintf('rawmodel_color_%s_size%d_ax%d.mat',folder,sxBySyTex,ax);
if ~exist(modelfile,'file')
    tic
    [template bx by] = learnColor(folder,sxBySyColor,ax,color_thres);
    disp(sprintf('selecting color features : %.3f seconds',toc));
else
    load(modelfile,'template','bx','by');
end

if isempty(template)
	colorSym = [];
else
    tx = bx; ty = by;
    Mx = template.Mx; My = template.My; gain = template.gain;

    colorSym = 0.8 * ones(bx,by,3) * 255;
    for iF = 1:length(gain)
        if sum(gain(iF)<color_thres)
            continue;
        end
        top = max( 1, Mx(iF) - ax );
        bottom = min( bx, Mx(iF) + ax );
        left = max( 1, My(iF) - ay );
        right = min( by, My(iF) + ay );
        vec = template.meanHist(iF,:);
        vec = vec;
        %vec = vec / sum(vec);
        colorSym(top:bottom,left:right,1) = vec(1);
        colorSym(top:bottom,left:right,2) = vec(2);
        colorSym(top:bottom,left:right,3) = vec(3);
    end
    
    %colorSym = template.mI;
    colorSym2 = template.varMap * 5e+4;
end





% ==========================================================
% illustrate the mixed template
% ==========================================================
basisSym = imresize(basisSym,[tx ty],'bilinear');
texSym = imresize(texSym,[tx ty],'bilinear');
flatSym = imresize(flatSym,[tx ty],'bilinear');
towrite = -basisSym;
towrite = 255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
imwrite(towrite,['sym_sk_' folder '.png']);
%imshow(basisSym);
towrite = -texSym;
towrite = 255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
imwrite(uint8(towrite),['sym_tex_' folder '.png']);
%imshow(texSym);
towrite = -flatSym;
towrite = 255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
if std(towrite) < 1e-3
    towrite = uint8( 255*ones(size(towrite)) );
end
imwrite(uint8(towrite)+128,['sym_flat_' folder '.png']);
%imshow(flatSym);
mixedSym = zeros(tx,ty,3);
mixedSym(:,:,2) = mixedSym(:,:,2) - basisSym ;
mixedSym(:,:,1) = mixedSym(:,:,1) - basisSym ;
mixedSym(:,:,3) = mixedSym(:,:,3) - basisSym ;
mixedSym(:,:,2) = mixedSym(:,:,2) - 0.8*texSym .* double(texSym>0);
mixedSym(:,:,3) = mixedSym(:,:,3) - 0.8*texSym .* double(texSym>0);
towrite = mixedSym;
towrite = 255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
imwrite(uint8(towrite),['sym_st _' folder '.png']);
mixedSym(:,:,3) = mixedSym(:,:,2) + 1*flatSym;
mixedSym(:,:,3) = mixedSym(:,:,1) + 1*flatSym;
mixedSym = 255 * (mixedSym-min(mixedSym(:)))/(max(mixedSym(:))-min(mixedSym(:)));
mixedSym = uint8(mixedSym);
imshow(mixedSym);
imwrite(mixedSym,['sym_' folder '.png']);
if ~isempty(colorSym)
	pause(1);
	imshow((colorSym));
    imwrite((colorSym),['symc_' folder '_ax' num2str(ax_color) '.png']);
    imwrite(uint8(min(255,colorSym2)),['varCI_' folder '_ax' num2str(ax_color) '.png']);
end


