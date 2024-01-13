clear all; close all; clc;

nTrials = 30;

s = rng;
v1 = normrnd(0,0.1,1,nTrials);
v2 = normrnd(0,0.1,1,nTrials);
w1 = zeros(1,nTrials);
w1(1) = 1;
w2 = zeros(1,nTrials);
w2(1) = 1;
for i = 2:nTrials
    w1(i) = w1(i-1) + v1(i-1);
    w2(i) = w2(i-1) + v2(i-1);
end
figure;
plot(0:nTrials-1,w1,'LineWidth',1.5)
hold on
plot(0:nTrials-1,w2,'LineWidth',1.5)
% ylim([0 2])
% xlim([0 20])
legend('W1','W2')
xlabel("trials")
ylabel("W")
title("Drift")