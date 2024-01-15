
%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_9


%% Paper DataSet
clear;clc;close all;

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\sparsenet');
rng shuffle
load('IMAGES.mat') 
load('IMAGES_RAW.mat')

A = rand(256) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));

figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(IMAGESr(:,:,i))
end
print -depsc Paper_W.eps
figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(IMAGES(:,:,i))
end
print -depsc Paper_unW.eps

figure(1), colormap(gray)
sparsenet

%% Caltech101 DataSet
clear;clc;close all;

Imagefolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\Caltech101';
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\Caltech101');
ImageNumbers=randi(9144,[1 10]);
files=dir(fullfile(Imagefolder,'*.jpg'));
[filenames{1:size(files,1)}] = deal(files.name);
for i=1:10
imgFile=filenames{1,ImageNumbers(i)};
ImageStruct{i}=im2double(rgb2gray(imread(imgFile)));
ImageSize(1,i)=size(ImageStruct{i},1);
ImageSize(2,i)=size(ImageStruct{i},2);
end
minImageSize=min(min(ImageSize));
for i=1:10
Image(:,:,i)=imresize(ImageStruct{i}, [minImageSize minImageSize]) ;
end
for i=1:10
X=Image(:,:,i);
mX = bsxfun(@minus,X,mean(X)); %remove mean
fX = fft(fft(mX,[],2),[],3); %fourier transform of the images
spectr = 0.25*sqrt(mean(abs(fX).^2)); %Mean spectrum
IMAGES(:,:,i) = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); 
end
save('Caltech101_unwhitened.mat','Image');
save('Caltech101_whitened.mat','IMAGES');

clear;clc;close all;
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\sparsenet');
rng shuffle
load('Caltech101_whitened.mat') 
load('Caltech101_unwhitened.mat')
A = rand(256) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));

figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(Image(:,:,i))
end
print -depsc Caltech101_unW.eps
figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(IMAGES(:,:,i))
end
print -depsc Caltech101_W.eps

figure(1), colormap(gray)
sparsenet

%% Yale DataSet
clear;clc;close all;

Imagefolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\yale';
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\yale');
ImageNumbers=randi(5400,[1 10]);
files=dir(fullfile(Imagefolder,'*.pgm'));
[filenames{1:size(files,1)}] = deal(files.name);
for i=1:10
imgFile=filenames{1,ImageNumbers(i)};
ImageStruct{i}=im2double(imread(imgFile));
ImageSize(1,i)=size(ImageStruct{i},1);
ImageSize(2,i)=size(ImageStruct{i},2);
end
minImageSize=min(min(ImageSize));
for i=1:10
Image(:,:,i)=imresize(ImageStruct{i}, [minImageSize minImageSize]) ;
end
for i=1:10
X=Image(:,:,i);
mX = bsxfun(@minus,X,mean(X)); %remove mean
fX = fft(fft(mX,[],2),[],3); %fourier transform of the images
spectr = 0.25*sqrt(mean(abs(fX).^2)); %Mean spectrum
IMAGES(:,:,i) = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); 
end
save('yale_unwhitened.mat','Image');
save('yale_whitened.mat','IMAGES');

clear;clc;close all;

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\sparsenet');
rng shuffle
load('yale_whitened.mat') 
load('yale_unwhitened.mat') 
A = rand(256) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));

figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(Image(:,:,i))
end
print -depsc Yale_unW.eps
figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(IMAGES(:,:,i))
end
print -depsc Yale_W.eps

figure(1), colormap(gray)
sparsenet


%% MNIST DataSet
clear;clc;close all;

Imagefolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\MNIST';
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\MNIST');
images = loadMNISTImages('t10k-images.idx3-ubyte');
ImageNumbers=randi(10000,[1 10]);
for i=1:10
ImageStruct{i}=reshape(images(1:784,ImageNumbers(i)),[28 28]);
ImageSize(1,i)=size(ImageStruct{i},1);
ImageSize(2,i)=size(ImageStruct{i},2);
end
minImageSize=min(min(ImageSize));
for i=1:10
Image(:,:,i)=imresize(ImageStruct{i}, [minImageSize minImageSize]) ;
end
for i=1:10
X=Image(:,:,i);
mX = bsxfun(@minus,X,mean(X)); %remove mean
fX = fft(fft(mX,[],2),[],3); %fourier transform of the images
spectr = 0.25*sqrt(mean(abs(fX).^2)); %Mean spectrum
IMAGES(:,:,i) = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); 
end
save('MNIST_unwhitened.mat','Image');
save('MNIST_whitened.mat','IMAGES');

clear;clc;close all;

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\sparsenet');
rng shuffle
load('MNIST_whitened.mat') 
load('MNIST_unwhitened.mat') 

A = rand(256) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));

figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(Image(:,:,i))
end
print -depsc MNIST_unW.eps
figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:10
subplottight(2,10,i)
imshow(IMAGES(:,:,i))
end
print -depsc MNIST_W.eps

figure(1), colormap(gray)
sparsenet

%% Part3_a
clear all; close all; clc;

video = read(VideoReader('BIRD.avi'));
temp=VideoReader('BIRD.avi')
frames = temp.NumberOfFrames;
for i = 1 : frames
    Video(:,:,i) = imcrop(im2double(rgb2gray(video(:, :, :, i))),[0,0,288,288]);
end
for j=1:12
    if j~=12
    FRAME=Video(:,:,1+(j-1)*10:10*j);
    else 
    FRAME=Video(:,:,1+(j-1)*10:10*j-2);
    end
    for i=1:size(FRAME,3)
      X=FRAME(:,:,i);
     mX = bsxfun(@minus,X,mean(X)); %remove mean
     fX = fft(fft(mX,[],2),[],3); %fourier transform of the images
 spectr = 0.35*sqrt(mean(abs(fX).^2)); %Mean spectrum
IMAGE(:,:,i) = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); 
    end
    IMAGEs.FRAME{j}=IMAGE;
end
save('Frames.mat','IMAGEs');

clear;clc;close all;

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\sparsenet');
rng shuffle
load('Frames.mat') 
IMAGES=IMAGEs.FRAME{1,2};

A = rand(64) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));

figure(1), colormap(gray)
sparsenet

%% Part3_b
clear;clc;close all;

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW09\Data\sparsenet');
rng shuffle
load('Frames.mat') 

A = rand(64) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));

patchSize = 10;
noise_var = 0.01;
beta = 2.2;
sigma = 0.316;
tol = 0.01;

for i=1:12
    IMAGES=IMAGEs.FRAME{1,i};
    X=[];
   for j=1:size(IMAGES,3)
       X=[X,patcher(IMAGES(:,:,j),8)];
   end
   S{i}=cgf_fitS(A,X,noise_var,beta,sigma,tol);
end

figure;
set(gcf,'units','points','position',[0,0,1600,650])
 patch = randi(12250);
    for i = 1:12
        subplot(3,4,i)
        plot(abs(S{i}(:,patch)));
        ylabel("Value")
        xlabel("Coefficient")
        if i~=12
        title(sprintf("Frame: %d to %d",1+(i-1)*10,10*i))
        else 
        title(sprintf("Frame: %d to %d",1+(i-1)*10,10*i-2))
        end    
    end
 print -depsc patch1.eps
   


function images = loadMNISTImages(filename)

fp = fopen(filename, 'rb');
assert(fp ~= -1, ['Could not open ', filename, '']);

magic = fread(fp, 1, 'int32', 0, 'ieee-be');
assert(magic == 2051, ['Bad magic number in ', filename, '']);

numImages = fread(fp, 1, 'int32', 0, 'ieee-be');
numRows = fread(fp, 1, 'int32', 0, 'ieee-be');
numCols = fread(fp, 1, 'int32', 0, 'ieee-be');

images = fread(fp, inf, 'unsigned char');
images = reshape(images, numCols, numRows, numImages);
images = permute(images,[2 1 3]);

fclose(fp);
% Reshape to #pixels x #examples
images = reshape(images, size(images, 1) * size(images, 2), size(images, 3));
% Convert to double and rescale to [0,1]
images = double(images) / 255;

end

function [X] = patcher(img,patchSize)
    sz = size(img);
    p = {};
    inds1 = 1:patchSize:sz(1)-patchSize;
    inds2 = 1:patchSize:sz(2)-patchSize;
    for i = 1:length(inds1)
        i1 = inds1(i);
        for j = 1:length(inds2)
            j1 = inds2(j);
            p{(i-1)*length(inds2) + j} = img(i1:i1+patchSize-1,j1:j1+patchSize-1);
            X(:,(i-1)*length(inds2) + j) = reshape(img(i1:i1+patchSize-1,j1:j1+patchSize-1),patchSize^2,1);
        end
    end
    
end
function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    ax = subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n])
    if(nargout > 0)
      h = ax;
    end
end

