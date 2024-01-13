clear
clc
close all
trial = 20;


r = ones(1,trial);
w0 = [0;0];
u = [ones(1,trial);ones(1,trial/2),zeros(1,trial/2)];
w_n = eye(2)*0.02;
tau = 1.2;
sigma0 = eye(2)*0.6;
sigma = zeros(2,2,trial);
sigma(:,:,1)=sigma0;
w0=[0;0];
w=zeros(2,20);
w(:, 1) = w0;
sigma1 = [];
sigma2 = [];
for i=1:trial-1
    sigma(:,:,i+1) = sigma(:,:,i)+ w_n;
    G = sigma(:,:,i+1)*u(:,i)*(u(:,i)'*sigma(:,:,i+1)*u(:,i)+tau^2)^(-1);
    sigma(:,:,i+1) =sigma(:,:,i+1)-G*u(:,i)'*sigma(:,:,i+1);
    w(:,i+1) = w(:,i)+G*(r(i)-u(:,i)'*w(:,i));

end

figure;
subplot(2,1,1)
plot(0:trial-1,w(1,:))
hold on
plot(0:trial-1,w(2,:))
xlim([0 trial-1])
hold on
title('Bachward Block                                                                                                                  Mean')
xlabel('trials');
ylabel('W');
legend('W1','W2');
subplot(2,1,2)
plot(1:trial,squeeze(sigma(1,1,:)))
hold on 
plot(11:trial,squeeze(sigma(2,2,11:trial)))  
xlim([1 trial-1])
ylim([0 1])
hold on
title('                                                                                                                      Variance')
xlabel('trials');
ylabel('W');
legend('sigma1','sigma2');



sigma1 = sigma(:,:,2);
sigma9 = sigma(:,:,12);
sigma19 = sigma(:,:,20);
sigmaVec = {sigma1, sigma9, sigma19};
w1 = w(:,2);
w9 = w(:,12);
w19 = w(:,20);
ww = {w1, w9, w19};
r = 0.3:0.3:3;
r = flip(r);
theta = 0:0.01:2*pi;
time = [2 12 20];
figure
for sub = 1:3
    subplot(1,3,sub)
    tmpW = ww{sub};
    sigmaTmp = sigmaVec{sub};
    for i = 1:length(r)
        x = r(i)*sin(theta);
        y = r(i)*cos(theta);
        for j = 1:length(x)
            tmp = [x(j);y(j)];
            tmp = sigmaTmp * tmp;
            x(j) = tmp(1)+tmpW(1);
            y(j) = tmp(2)+tmpW(2);
        end
        scatter(x,y,'k')
        plot(x,y,'r')
        axis square
        hold on
    end
    mu = ww{sub}';
    sigma = sigmaVec{sub};
    hold on
    xlim([-1 2])
    ylim([-1 2])
    xlabel('W1')
    ylabel('W2')
    title(sprintf("t = %d",time(sub)))
end
