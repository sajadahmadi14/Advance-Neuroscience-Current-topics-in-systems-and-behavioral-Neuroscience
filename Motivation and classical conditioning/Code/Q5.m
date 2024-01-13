clear
clc
close all
nTrials = 22;


r = [ones(1,nTrials/2-2),-1*ones(1,nTrials/2+2)];
w0 = [0;0];
u = ones(1,nTrials);
w_n = 0.08;
tau = 0.5;

sigma = zeros(1,1,nTrials);
sigma0 = 0.6;
sigma(:,:,1)=sigma0;
w0=[0;0];
w=zeros(2,20);
w(:, 1) = w0;

for i=1:nTrials-1
    sigma(:,:,i+1) = sigma(:,:,i)+ w_n;
    G = sigma(:,:,i+1)*u(:,i)*(u(:,i)'*sigma(:,:,i+1)*u(:,i)+tau^2)^(-1);
    sigma(:,:,i+1) =sigma(:,:,i+1)-G*u(:,i)'*sigma(:,:,i+1);
    w(:,i+1) = w(:,i)+G*(r(i)-u(:,i)'*w(:,i));

end


figure;
subplot(2,1,1)
plot(0:nTrials-1,w(1,:))


xlim([0 20])
ylim([-1.2 1.2])
title("s1 -> r      s1 -> -r     mean")
ylabel('W')
xlabel('Trials')

subplot(2,1,2)
plot(0:nTrials-1,squeeze(sigma(1,1,:)))


xlim([0 20])
ylim([0 1])
title(' Variance')
ylabel('sigma')
xlabel('Trials')

