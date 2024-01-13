% clear all; close all; clc;
addpath('D:\Eye tracking database\matlabPyrTools-master\matlabPyrTools-master');
addpath('D:\Eye tracking database\SaliencyToolbox2.3\SaliencyToolbox')
addpath('D:\Eye tracking database\voc-release-3.1-win-master\voc-release-3.1-win-master')
addpath('D:\Eye tracking database\FaceDetect')
addpath('D:\Eye tracking database\LabelMeToolbox-master')
addpath('D:\Eye tracking database\JuddSaliencyModel\JuddSaliencyModel')
addpath('D:\Eye tracking database\JuddSaliencyModel\JuddSaliencyModel\horizon code')
addpath('D:\Eye tracking database\LabelMeToolbox-master\LabelMeToolbox-master\features')
addpath('D:\Eye tracking database\LabelMeToolbox-master\LabelMeToolbox-master\imagemanipulation')
addpath('D:\Eye tracking database\JuddSaliencyModel')
addpath('D:\Eye tracking database\JuddSaliencyModel\JuddSaliencyModel\FelzenszwalbDetectors')
addpath('D:\Eye tracking database\Eye tracking database\Eye tracking database\ALLSTIMULI')
addpath('D:\Eye tracking database\FaceDetect')
addpath('D:\Eye tracking database\matlabPyrTools-master\matlabPyrTools-master\MEX')
datafolder = 'D:\Eye tracking database\Eye tracking database\Eye tracking database\DATA\krl';
stimfolder = 'D:\Eye tracking database\Eye tracking database\Eye tracking\database\ALLSTIMULI';
files=dir(fullfile(datafolder,'*.mat'));
[filenames{1:size(files,1)}] = deal(files.name);
Nstimuli = size(filenames,2);


numImg = 30;
for i=21:numImg
    load(fullfile(datafolder,filenames{i}))
    stimFile = eval([filenames{i}(1:end-4)]);
    imgName = stimFile.imgName;
    map = {};
    load model
    img = imread(imgName);
       k(:,:,i)=salineccy(img);
% map = {};
% for i = 1:12
%     a=IMAGES_holder{i};
%     
%    for k=1:10
%     img = a(:,:,k);
%     img = repmat(img,1,1,3);
%     saliencyMap(:,:,k)=salineccy(img);
%    end
%     map{i} = saliencyMap;
end
save map.mat k
%% Patching
clear all; close all; clc;
load map.mat
load IMAGES.mat
pSize = 8;
p = [];
th = (pSize*pSize)/2;
th = 30;
cnt = 1;
b=1;
n=1;
k=zeros(128,128);
for i = 1:10
    n=1;
    img = IMAGES(:,:,i);
    salMap = map{i};
    step = pSize - pSize/2;
    for r = 1:step:size(img,1)-step
        for c = 1:step:size(img,2)-step
           % pTmp = img(r:r+pSize-1,c:c+pSize-1);
            %pTmp = reshape(pTmp,pSize^2,1);
            sc = sum(salMap(r:r+pSize-1,c:c+pSize-1),'all');
            if (sc >= th)
             %   p(:,cnt) = pTmp;
             k(b,n,i)=1;   
            end
            cnt = cnt + 1;
            b=b+1;
        end
        n=n+1;
        b=1;
    end
%     im{i}=p;
%     cnt=1;
%     p=[];
end
save k k


 addpath('C:\Users\Pars\Desktop\sparsenet\sparsenet')
 addpath('C:\Users\Pars\Desktop\codee')
load IMAGES.mat
load imgr.mat
for s=1:size(k,3)
   crop(:,:,s) = imcrop(k(:,:,s),[0,0,768,768]);
end
IMAGES=crop;
A = rand(256) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));
figure(1), colormap(gray)
sparsenet



%% Functions
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

function IMAGES = whitening(images)
    N = size(images,1);
    imgNum = size(images,3);
    M = imgNum;
    [fx, fy]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
    rho=sqrt(fx.*fx+fy.*fy);
    f_0=0.4*N;
    filt=rho.*exp(-(rho/f_0).^4);
    for i=1:M
        image=images(:,:,i);
        If=fft2(image);
        imagew=real(ifft2(If.*fftshift(filt)));
        IMAGES(:,:,i) = imagew;
    end
    IMAGES=sqrt(0.1)*IMAGES./sqrt(mean(var(IMAGES)));
end


function saliencyMap=salineccy(img)
% origimgsize = size(img);
% imgName = stimFile.imgName;
% img = imread(imgName);
[w, h, c] = size(img);
dims = [200, 200];
    FEATURES(:, 1:13) = findSubbandFeatures(img, dims);
    FEATURES(:, 14:16) = findIttiFeatures(img, dims);
    FEATURES(:, 17:27) = findColorFeatures(img, dims);
    FEATURES(:, 28) = findTorralbaSaliency(img, dims);
    %FEATURES(:, 29) = findHorizonFeatures(img, dims);
    %FEATURES(:, 30:31) = findObjectFeatures(img, dims);
    %FEATURES(:, 32) = findDistToCenterFeatures(img, dims);
lengthF = length((FEATURES(1,:)));
load model
    meanVec = model.whiteningParams(1, lengthF);
    stdVec = model.whiteningParams(2, lengthF);
    FEATURES=FEATURES-repmat(meanVec, [size(FEATURES, 1), 1]);
    FEATURES=FEATURES./repmat(stdVec, [size(FEATURES, 1), 1]);
    saliencyMap = (FEATURES*model.w(1:lengthF)') + model.w(end);
    saliencyMap = (saliencyMap-min(saliencyMap))/(max(saliencyMap)-min(saliencyMap));
    saliencyMap = reshape(saliencyMap, dims);
    saliencyMap = imresize(saliencyMap, [w, h]);
  

end