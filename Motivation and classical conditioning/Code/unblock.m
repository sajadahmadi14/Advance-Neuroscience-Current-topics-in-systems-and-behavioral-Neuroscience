clear
clc
close all
nTrials = 30;


r = [ones(1,nTrials/2),2*ones(1,nTrials/2)];
w0 = [0;0];
u = [ones(1,nTrials);zeros(1,nTrials/2),ones(1,nTrials/2)];
w_n = eye(2)*0.01;
tau = 0.6;
sigma0 = eye(2)*0.6;
%[w,sigma] = myKalman(W_noise,w0,nTrials,tau,sigma0,C,r);
sigma = zeros(2,2,nTrials);
sigma0 =0.6 * eye(2);
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
hold on
plot(0:nTrials-1,w(2,:))
hold on
title('UnBlocking                                                                                                                  Mean')
xlabel('trials');
ylabel('W');
legend('W1','W2');
subplot(2,1,2)
plot(1:nTrials,squeeze(sigma(1,1,:)))
hold on 
plot(15:nTrials,squeeze(sigma(2,2,15:nTrials)))  
hold on
title('                                                                                                                      Variance')
xlabel('trials');
ylabel('W');
legend('sigma1','sigma2');
