%% Question 6:
clear; clc; close all;

thrU = 1;
thrD = -1;
x0 = 0;
dt = 0.01;
sigma = 1;
Bias = 0.1;
mean = 0;

figure;
for k = 1:10
    [choice(k), RT(k), x{k}] = two_choice_trial(x0, Bias, dt, mean, sigma, thrU, thrD);
    if choice(k)==1
    plot([0: dt: RT(k)-dt], x{k},'r');
    hold on;
    else
         plot([0: dt: RT(k)-dt], x{k},'k');
    end
end
title('two choice 1 -1')
xlabel('time')
% hold all;
ylabel('X for 10 experiment')
% hold on; 
% plot([0: dt: RT(k)-dt], x{k})
axis tight
function [choice, RT, x] = two_choice_trial(x0, B, dt, m, sigma, thrU, thrD)
    x(1) = x0;
    i = 2;
    while 1
        x(i) = x(i-1) + B*dt + sigma * normrnd(m, sqrt(dt), 1);
        if x(i) > thrU
            RT = i * dt;
            choice = 1;
            x(i) = thrU;
            break;
        end
        if x(i) < thrD
            RT = i * dt;
            choice = -1;
            x(i) = thrD;
            break;
        end
        i = i + 1;
    end
end

