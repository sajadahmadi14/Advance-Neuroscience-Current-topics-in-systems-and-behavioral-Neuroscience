clear all; close all; clc;

trlNum = 100;
posTh = 5;
negTh = -5;
bias = 0.1;
sigma = 1;
x0 = 0;
dt = 0.01;
RTVec = [];
choiceVec = [];
for trl = 1:trlNum
    trl
    [RT,choice] = two_choice_trial(posTh,negTh,sigma,x0,bias,dt);
    RTVec(trl) = RT;
    choiceVec(trl) = choice;
end

crctInd = find(choiceVec == 1);
wrInd = find(choiceVec == -1);
nbin = 10;
[cnt,x] = hist(RTVec,nbin);
[cntCorrect,xCorrect] = hist(RTVec(crctInd),nbin);
[cntWR,xWR] = hist(RTVec(wrInd),nbin);

figure;
subplot(1,3,1)
plot(x,cnt,'k','LineWidth',1.5)
xlim([0,max(RTVec)])
ylim([0,700])
xlabel("RT(s)")
ylabel("Count")
title("RT distribution")

subplot(1,3,2)
plot(xCorrect,cntCorrect,'k','LineWidth',1.5)
xlim([0,max(RTVec)])
ylim([0,700])
xlabel("RT(s)")
ylabel("Count")
title("RT distribution for correct trials")

subplot(1,3,3)
plot(xWR,cntWR,'k','LineWidth',1.5)
xlim([0,max(RTVec)])
ylim([0,700])
xlabel("RT(s)")
ylabel("Count")
title("RT distribution for  wrong trials")

figure;
cntCorrect = cntCorrect/max(cntCorrect);
cntWR = cntWR/max(cntWR);
plot(xCorrect,cntCorrect,'k','LineWidth',1.5)
hold on
plot(xWR,cntWR,'r','LineWidth',1.5)
xlabel("RT(s)")
ylabel("Count")
title("Normalized Count")
legend("Correct","Wrong")
xlim([0,max(RTVec)])

function [RT,choice] = two_choice_trial(posTh,negTh,sigma,x0,bias,dt)
    x = x0;
    t = 0;
    while (x <= posTh && x >= negTh)
        x = x + bias * dt + sigma * normrnd(0,sqrt(dt));
        t = t+dt;
    end
    RT = t;
    if (x >= posTh)
        choice = 1;
    else
        choice = -1;
    end
end