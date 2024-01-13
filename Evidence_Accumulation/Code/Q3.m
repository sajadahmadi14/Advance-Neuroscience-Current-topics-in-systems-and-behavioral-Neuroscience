clear;
clc;
close all;
trials = 1000;
sigma = 1;
Bias = 0.1;
mean = 0;
dt = 0.1;
x0 = 0;
TVec = linspace(0.5, 10, 100);

for k = 1: length(TVec)
    for itrial = 1: trials
        [choice(k, itrial), x] = sm(Bias, sigma, mean, dt, x0, TVec(k));
    end
end

for k = 1: length(TVec)
    errorP(k) = length(find(choice(k, :) == -1))/length(choice(k, :));
end

clear sigma B
pth_fun = @(t, sigma, T, B) (1/sqrt(2*pi*sigma*T))*exp(-((t - B*T).^2)/(2*sigma*T));
for k = 1: length(TVec)
    errorTh(k) = integral(@(t)pth_fun(t, 1, TVec(k), 0.1), -inf, 0);
end
figure; 
scatter(TVec, errorP); 
hold on;
plot(TVec, errorTh);
title('Error rate'); xlabel('Time'); ylabel('Error');
legend('Simulation', 'Theory')
axis tight


function [choice, x] = sm(B, sigma, mean, dt, x0, T)
    x(1) = x0;
    t=T/dt;
    for i = 2:t
        x(i) = x(i-1) + B*dt + sigma * normrnd(mean, sqrt(dt)); 
    end
    choice = sign(x(end));
end

 