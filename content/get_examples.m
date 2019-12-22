function I = get_examples(folder)

imageNames = dir([folder '/*.jpg']);
if isempty(imageNames)
    imageNames = dir([folder '/*.bmp']);
end
if isempty(imageNames)
    imageNames = dir([folder '/*.png']);
end
nimage = length(imageNames);
I = cell(1,nimage);
for i = 1 : nimage
    im = imread([folder '/' imageNames(i).name]);
    if size(im,3)==3
        im = rgb2gray(im);
    end
    I{i} = double(im);
end
