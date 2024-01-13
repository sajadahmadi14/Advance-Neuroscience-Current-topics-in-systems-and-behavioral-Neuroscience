
clear all; close all; clc;
rng shuffle

Bias = 0.1;
sigma = 1;
dt = 0.1;
T = 10;
trials = 1000;
xVec = [];
for trl = 1:trials
    [~,x] = sm(Bias,sigma,0,dt,0,T);
    xVec = [xVec;x];
end
len = round(T/dt);
t = linspace(0,T,len);
stdX = std(xVec,[],1);
meanX = mean(xVec,1);
meanSim = Bias*t;
stdSim = sqrt(sigma*t);

figure;
plot(t,xVec(1,:));
hold on
f=rand(1,50);
f=f*100;
f=round(f);
for j=f
plot(t,xVec(j,:));
hold on
end
plot(t,meanX,'k','linewidth',1.5)
plot(t,meanX+stdX,'--k','linewidth',1.5)
plot(t,meanX-stdX,'--k','linewidth',1.5)
plot(t,meanSim,'r')
plot(t,meanSim+stdSim,'--r')
plot(t,meanSim-stdSim,'--r')

title("Mean and STD of X")
xlabel("Time")
ylabel("X")
legend("x","Simulation Mean and std","Theoretical Mean and std")


function [choice, x] = sm(B, sigma, mean, dt, x0, T)
    x(1) = x0;
    t=T/dt;
    for i = 2:t
        x(i) = x(i-1) + B*dt + sigma * normrnd(mean, sqrt(dt)); 
    end
    choice = sign(x(end));
end

 