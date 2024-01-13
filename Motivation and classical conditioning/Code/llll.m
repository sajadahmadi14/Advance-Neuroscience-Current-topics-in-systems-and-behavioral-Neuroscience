clear all; close all; clc;

nTrials = 100;
s = rng;
v = normrnd(0,0.1,1,nTrials);
vr = normrnd(0,0.5,1,nTrials);
phi = normrnd(0,2,1,nTrials);

c = zeros(1,nTrials);
c(40) = 1;
c(90) = 1;
phi(40) = -2;
phi(90) = 4;
ind = find(c==1);
w = zeros(1,nTrials);
r = w;
w(1) = 0;
r(1) = w(1) + vr(1) + c(1)*phi(1);
for i = 2:nTrials
    w(i) = w(i-1) + v(i-1) + c(i-1) * phi(i-1);
    r(i) = w(i) + vr(i) + c(i)*phi(i);
end
figure;
subplot(3,1,1)
scatter(0:nTrials-1,w,20,'filled','k')
hold on
scatter(0:nTrials-1,r,'xk')
hold on
scatter([0,ind],[w(1),w(ind+1)],500,'k','LineWidth',2)
hold on
for i = 1:length(ind)
    xline(ind(i),':k');
end
xlim([0 nTrials])
ylim([-4 4])
title("Slow Draft and Dramatic Changes")
xlabel("t")
ylabel("w")
legend("w(t)","r(t)",'location','northwest')

% Learning Weights
gamma = 3.3;
W_noise = 0.01;
w0 = 0;
tau = 0.7;
sigma0 = 0.6;
u = ones(1,nTrials);
    w = zeros(1,nTrials);
    sigma = zeros(1,nTrials);
    betaVec = zeros(1,nTrials);
    w(1) = w0;
    sigma(1) = sigma0;
    for i = 2:nTrials
        % prediction
        sigmap = sigma(i-1) + W_noise;
        
        % update
        G = sigmap * u(i) * (u(i)*sigmap + tau^2)^-1;
        sigma(i) = sigmap - G*u(i)*sigmap;
        w(i) = w(i-1) + G*(r(i) - u(i)*w(i-1));
        
        % NE
        beta = (r(i) - u(i)*w(i))^2 / (u(i)*sigma(i) + tau^2);
        betaVec(i) = beta; 
        if beta > gamma
            sigma(i) = 100;
        end
    end
subplot(3,1,2)
scatter(0:nTrials-1,w,20,'filled','k')
hold on
scatter(0:nTrials-1,r,'xk')
hold on
scatter(0:nTrials-1,w,'k');
xlim([0 nTrials])
ylim([-4 4])
legend("theoretical w(t)","r(t)","estimated w(t)",'location','northwest')
title(sprintf("Learned w alongside r and theoretical w, gamma = %.2f",gamma))
xlabel("t")

subplot(3,1,3)
plot(sigma,'--k','LineWidth',1.5)
hold on
plot(betaVec,'k','LineWidth',1.5)
hold on
yline(gamma,'-.k','LineWidth',1.5);
ylim([0 10])
legend("ACh","NE","$gamma$",'Interpreter','LaTex','location','northwest')
xlabel("t")

n = 100;
gammaVec = linspace(0,30,n);
mse = zeros(1,n);
for j = 1:1000
    nTrials = 100;
    s = rng;
    v = normrnd(0,0.1,1,nTrials);
    vr = normrnd(0,0.5,1,nTrials);
    phi = normrnd(0,2,1,nTrials);
    
    c = zeros(1,nTrials);
    c(40) = 1;
    c(90) = 1;
    phi(40) = -2;
    phi(90) = 4;
    ind = find(c==1);
    w = zeros(1,nTrials);
    r = w;
    w(1) = 0;
    r(1) = w(1) + vr(1) + c(1)*phi(1);
    for i = 2:nTrials
        w(i)=w(i-1)+v(i-1)+c(i-1)*phi(i-1);
        r(i)=w(i)+vr(i)+c(i)*phi(i);
    end
    mseTmp = [];
    for i1 = 1:n
        gamma = gammaVec(i1);
           w1 = zeros(1,nTrials);
    sigma = zeros(1,nTrials);
    betaVec = zeros(1,nTrials);
    w1(1) = w0;
    sigma(1) = sigma0;
    for i = 2:nTrials
        % prediction
        sigmap = sigma(i-1) + W_noise;
        
        % update
        G = sigmap * u(i) * (u(i)*sigmap + tau^2)^-1;
        sigma(i) = sigmap - G*u(i)*sigmap;
        w1(i) = w1(i-1) + G*(r(i) - u(i)*w1(i-1));
        
        % NE
        beta = (r(i) - u(i)*w1(i))^2 / (u(i)*sigma(i) + tau^2);
        betaVec(i) = beta; 
        if beta > gamma
            sigma(i) = 100;
        end
    end
        tmp = sum((w1-w).^2);
        mseTmp = [mseTmp,tmp];
    end
    mse = mse + mseTmp;
end
mse = mse/1000;
figure
plot(gammaVec,mse,'k','LineWidth',1.5)
title("MSE for different gamma")
xlabel("$\gamma$",'interpreter','LaTex')
ylabel("MSE")

