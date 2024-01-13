clear; clc; close all;

x01 = 0.5;
x02 = 0;
Bias1 = 0.05;
Bias2 = 0.05;
sigma1 = 1;
sigma2 = 1;
thr1 = 2; 
thr2 = 2;
dt = 0.1;
mean = 0;

figure; 
for k = 1: 8
    
    [choice(k),RT(k), x1{k}, x2{k}, winner(k)] = race_trial(x01, x02,Bias1, Bias2, sigma1, sigma2, thr1, thr2, mean, dt);
    subplot(2,4,k)
    plot([0: dt: RT(k)-dt], x1{k},'k', 'linewidth', 1.5);
    hold all;
    plot([0: dt: RT(k)-dt], x2{k}, 'r', 'linewidth', 1.5);
    hold all;
    yline(thr2, 'k', 'linewidth', 4);
    xlabel('Time'); ylabel('x(t)');
    title('Race X1 Vs X2')
    legend('x1','x2')
end



function [choice, RT, x1, x2, winner] = race_trial(x01, x02, B1, B2, sigma1, ...
    sigma2, thr1, thr2, m, dt)
    x1(1) = x01;
    x2(1) = x02;
    i = 2;
    while 1
        x1(i) = x1(i-1) + B1*dt + sigma1 * normrnd(m, sqrt(dt), 1);
        x2(i) = x2(i-1) + B2*dt + sigma2 * normrnd(m, sqrt(dt), 1);
        if x1(i) > thr1
            RT = i * dt;
            choice = 1;
            x1(i) = thr1;
            winner = 1;
            break;
        end
        if x2(i) > thr2
            RT = i * dt;
            choice = 1;
            x2(i) = thr2;
            winner = 2;
            break;
        end
        i = i + 1;
    end
end

