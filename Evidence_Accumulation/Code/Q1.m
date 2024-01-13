clear; clc; close all;
sigma = 1;
Bias = 1;
mean = 0;
dt = 0.01;
x0 = 0;
T = 10;

[choice, x] = sm(Bias, sigma, mean, dt, x0, T);

figure; plot([0: dt: T-dt], x);
xlabel('Time'); ylabel('x');
title('Drift diffusion model');


function [choice, x] = sm(B, sigma, mean, dt, x0, T)
    x(1) = x0;
    t=T/dt;
    for i = 2:t
        x(i) = x(i-1) + B*dt + sigma * normrnd(mean, sqrt(dt)); 
    end
    choice = sign(x(:,T/dt));
end

 