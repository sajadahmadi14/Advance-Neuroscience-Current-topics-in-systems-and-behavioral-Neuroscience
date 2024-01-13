
clear all; close all; clc;
addpath('C:\Users\Sajad\Desktop\Data')
obj = VideoReader('BIRD.avi');
vid = read(obj);
frames = obj.NumberOfFrames;
for i = 1 : frames
    a = rgb2gray(vid(:, :, :, i));
    a = im2double(a);
    crop = imcrop(a,[0,0,280,280]);
    Vid(:,:,i) = crop;
end
save Vid 

IMAGES_holder = {};
for i = 1:10:110
    a = Vid(:,:,i:i+9);
    images=a;
    N = size(images,1);
    imgNum = size(images,3);
    M = imgNum;
    [fx, fy]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
    rho=sqrt(fx.*fx+fy.*fy);
    f_0=0.4*N;
    filt=rho.*exp(-(rho/f_0).^4);
    for g=1:M
        image=images(:,:,g);
        If=fft2(image);
        imagew=real(ifft2(If.*fftshift(filt)));
        IMAGES(:,:,g) = imagew;
    end
    IMAGES=sqrt(0.1)*IMAGES./sqrt(mean(var(IMAGES)));
    a=IMAGES;
    IMAGES_holder{round(i/10+1)} = a;
end

save IMAGES_holder IMAGES_holder

A=rand(64)-0.5;
A=A*diag(1./sqrt(sum(A.*A)));
load IMAGES_holder.mat

patchSize=8;

for f = 1:10
    images = IMAGES_holder{f};
    X = [];
    
        for j = 1:10
            img = images(:,:,j);
            X = [X,patcher(img,patchSize)];
        end
        S{f} = cgf_fitS(A,X);
end
load s
for j = 1:4
    subplot(2,2,j)
    x = [];
    p = randi(9240);
    for i = 1:10
        a = S{i};
        x = [x,abs(a(:,p))];
        m=max(abs(a(:,p)));
        plot(abs(a(:,p))+(i),'k');
        yline((i),'b','linewidth',0.5);
        hold on
    end
    ylim tight
    xlim tight

    title(sprintf("Patch %d",p))
    ylabel("Time")
    xlabel("Coefficient")
%     axis tight
end



x=[];
figure;
for i = 1:10
    tmp = S{i};
    avg = mean(abs(tmp),2);
    x = [x,avg];
    plot(avg+i,'k');
    yline((i),'b');
    hold on
    axis tight
end
plot(mean(x,2)+(i+1),'k','linewidth',1.5)
   title('Average on all patches')
    ylabel("Time")
    xlabel("Coefficient")
Histograms
figure;
for i = 6:10
    a = S{i};
    subplot(2,5,i);
   histogram(a)  
   title(sprintf(" %d",i))
   sgtitle('coefficients')
    subplot(2,5,i-5);
    plot(a)  
axis tight
end
% figure;

% for i = 6:10
%     a = S{i};
%     subplot(5,1,i);
%    plot(a,'c')  
% axis tight
%     title(sprintf("Frame set %d",i))
% end
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

