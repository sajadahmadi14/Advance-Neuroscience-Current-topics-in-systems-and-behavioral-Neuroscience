
clear; clc; close all;
Bias = 0.1;
T = 10;
sigma = 1;
startV= -10:0.1:10;
for k = 1: length(startV)
    for p = 1: 1000
        choice(k, p) = sm2(startV(k), Bias, T, sigma);
    end
    p1(k) = length(find(choice(k, :) == -1))/length(choice(k, :));
end

figure;
cla()
k=scatter(startV, p1);

hold on
plot(startV, p1)
% Plot line
k2=lsline 
xlabel('Starting Point'); ylabel('Error rate');
title('simple model 2')
function [choice] = sm2(start_point, B, T, sigma)
    p = cdf('Normal', start_point, start_point + B*T, sqrt(sigma*T));
    p_decide = rand(1);
    if p_decide < p, choice = -1; else, choice = 1; end
end



 
 