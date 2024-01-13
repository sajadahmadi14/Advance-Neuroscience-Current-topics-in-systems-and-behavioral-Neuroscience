% clc
% clear all
% close all

addpath('C:\Users\Pars\Desktop\Dataa')
addpath('C:\Users\Pars\Desktop')
% load mnist
% load IMAGES_RAW
ff=im2double(ff);
IMAGESS=(ff(:,472,1));
N = size(IMAGESS,1);

% imgNum = size(C,1);
M = size(IMAGESS,3);



[fx, fy]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
rho=sqrt(fx.*fx+fy.*fy);
f_0=0.4*N;
filt=rho.*exp(-(rho/f_0).^4);

for i=1:M
    image=IMAGESS(:,:,i);
    If=fft2(image);
    imagew=real(ifft2(If.*fftshift(filt)));
%     IMAGES(:,i)=reshape(imagew,N^2,1);
    IMAGES(:,:,i) = imagew;
end
IMAGES=sqrt(0.1)*IMAGES./sqrt(mean(var(IMAGES)));
% save mnistWH IMAGES
% save mnistorg IMAGESS

