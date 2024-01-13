clc 
clear all
close all
load IMAGES_RAW.mat
load k 

b=1;
n=1;
pSize = 8;
step=4;
for o=1:10
 subplot(2,5,o)
 img=k(:,:,o); imshow(img)
 sgtitle('Whitened Orginal Images')
end
figure
for o=1:10
img=IMAGES(:,:,o);
n=1;
 for r = 1:step:size(img,1)-step
        for c = 1:step:size(img,2)-step
            if k(b,n,o)==0
           img(r:r+pSize-1,c:c+pSize-1)=0;
            end
            b=b+1;
        end
        n=n+1;
        b=1;
 end
 subplot(2,5,o)
 imshow(img)
 IMAGES(:,:,o)=img;
 sgtitle('Whitened Salient Parts TH=20')
 
end
 save imgr IMAGES