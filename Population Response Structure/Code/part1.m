
%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_1

%% Part1-section(a)

clear;
clc;
close all;

r=100; %% frequency 
SimuTime=1; 
deltaTau = 1/1000; 
NumTrials=150;
t = 0:deltaTau:SimuTime-deltaTau;
Spike=SpikePoissonGen(r,SimuTime,deltaTau,NumTrials);
figure;
RasterPlot(Spike,t);
xlabel('Time (s)','interpreter', 'latex');
ylabel('Trials','interpreter', 'latex');
title(['Raster Plot'],'fontsize',8,'interpreter', 'latex')

print -depsc fig1.eps


%% Part1-section(b)
r=100; %% frequency 
SimuTime=1; 
deltaTau = 1/1000; 
NumTrials=5000;
t = 0:deltaTau:SimuTime-deltaTau;
Spike=SpikePoissonGen(r,SimuTime,deltaTau,NumTrials);
SpikeCount=sum(Spike,2);
x=min(SpikeCount):1:max(SpikeCount);
h = hist(SpikeCount , length(x));
figure;
bar(x,h/NumTrials,'b');
hold on
x=min(SpikeCount):1:max(SpikeCount);
PoissonTheory=poisspdf(x,r);
plot(x,PoissonTheory,'r','linewidth',2);
xlabel('Count','interpreter', 'latex');
ylabel('Probability','interpreter', 'latex');
title(['Spike Count Probability Histogram'],'fontsize',8,'interpreter', 'latex')
text(115,0.02,'$$Poisson(\lambda=r=100)$$','interpreter', 'latex','color','r')

print -depsc fig2.eps

%% Part1-section(c)

[ISI,SpikesTime]=myISI(Spike);

x=min(ISI):0.001:max(ISI);
h = hist(ISI , length(x));
figure;
bar(x,h/(NumTrials*r),'b');
hold on
x=0:0.001:max(ISI);
expTheory=exppdf(x,1/r)/1000;
plot(x,expTheory,'r','linewidth',2);
xlabel('ISI (s)','interpreter', 'latex');
ylabel('Probability','interpreter', 'latex');
title(['Inter-Spike Interval(ISI) Histogram'],'fontsize',8,'interpreter', 'latex')
text(0.02,0.03,'$$exponential(\lambda=r=100)$$','interpreter', 'latex','color','r')
print -depsc fig3.eps

CV=std(ISI)/mean(ISI);


%% Part(a,b,c) for Renewal Process

% Part(a)
%k=randi([1 10],1);
k=5;
[RenewalSpike]=myRenewal(Spike,k);

% Part(b)
RenewSpikeCount=sum(RenewalSpike,2);
x=min(RenewSpikeCount):1:max(RenewSpikeCount);
h = hist(RenewSpikeCount , length(x));
figure;
bar(x,h/NumTrials);
hold on
xlabel('Count','interpreter', 'latex');
ylabel('Probability','interpreter', 'latex');
title(['ISI Histogram'],'fontsize',8,'interpreter', 'latex')
text(0.08,0.01,'$$Gamma(k=5,\lambda=r=100)$$','interpreter', 'latex','color','r')

print -depsc fig4.eps

% Part(c)
[RenewISI,SpikesTime]=myISI(RenewalSpike);

x=min(RenewISI):0.001:max(RenewISI);
h = hist(RenewISI , length(x));
figure;
bar(x,h/(NumTrials*r/k));
hold on
GammaTheory=gampdf(x,k,1/r)/1000;
plot(x,GammaTheory,'r','linewidth',2);
xlabel('ISI (s)','interpreter', 'latex');
ylabel('Probability','interpreter', 'latex');
title(['ISI Histogram'],'fontsize',8,'interpreter', 'latex');
text(0.08,0.01,'$$Gamma(k=5,\lambda=r=100)$$','interpreter', 'latex','color','r')

print -depsc fig5.eps

RenewCV=std(RenewISI)/mean(RenewISI);

%% Part1-section(g)
clear;
clc;
close all;

figure;
Deltat=0.001:10^-4:0.030;
Nth=[1 4 51];
t0=0.001;
for i=1:length(Nth)
Cv=(1/sqrt(Nth(i)))*((Deltat-t0)./Deltat);
scatter(Deltat,Cv)
hold on
plot(Deltat,1./sqrt(Nth(i))*ones(1,length(Deltat)),'color','r')
hold on
end
ylim([0 1.1])
xlabel('$$\overline{\Delta t}(s)$$','interpreter', 'latex');
ylabel('CV','interpreter', 'latex');
title(['figure 6 (softky \& koch 1993)'],'fontsize',8,'interpreter', 'latex');
text(0.015,1.04,'$$N_{th}=1 , t_0=0$$','interpreter', 'latex');
text(0.025,0.9,'$$N_{th}=1 , t_0=1ms$$','interpreter', 'latex');
text(0.015,0.53,'$$N_{th}=4 , t_0=0$$','interpreter', 'latex');
text(0.025,0.42,'$$N_{th}=4 , t_0=1ms$$','interpreter', 'latex');
text(0.015,0.17,'$$N_{th}=51 , t_0=0$$','interpreter', 'latex');
text(0.025,0.09,'$$N_{th}=51 , t_0=1ms$$','interpreter', 'latex');
print -depsc fig7.eps

% Simulation
Deltat=[0.001:0.0003:0.030];
FiringRate=1./Deltat;
Nth=[1 4 51];
t0=0.001;
deltaTau=0.001;
SimuTime=10;
NumTrials=2000;
figure;
for i=1:3
    for j=1:length(FiringRate)
        [Spike,t]=SpikePoissonGen(FiringRate(j),SimuTime,deltaTau,NumTrials);
        [RenewalSpike]=myRenewal(Spike,Nth(i));
        [RenewISI]=myISI(RenewalSpike);
        cv(i,j)=std(RenewISI+t0)/mean(RenewISI+t0)
    end
    scatter(Deltat,cv(i,:))
    hold on
    plot(Deltat,1./sqrt(Nth(i))*ones(1,length(Deltat)),'color','r')
    hold on
end
ylim([0 1.1])
xlabel('$$\overline{\Delta t}(s)$$','interpreter', 'latex');
ylabel('CV','interpreter', 'latex');
title(['figure 6 (softky \& koch 1993)'],'fontsize',8,'interpreter', 'latex');
text(0.015,1.04,'$$N_{th}=1 , t_0=0$$','interpreter', 'latex');
text(0.025,0.9,'$$N_{th}=1 , t_0=1ms$$','interpreter', 'latex');
text(0.015,0.53,'$$N_{th}=4 , t_0=0$$','interpreter', 'latex');
text(0.025,0.42,'$$N_{th}=4 , t_0=1ms$$','interpreter', 'latex');
text(0.015,0.17,'$$N_{th}=51 , t_0=0$$','interpreter', 'latex');
text(0.025,0.09,'$$N_{th}=51 , t_0=1ms$$','interpreter', 'latex');

print -depsc fig8.eps


%% functions
function [Spike,t]=SpikePoissonGen(r,SimuTime,deltaTau,NumTrials)
NumBins=floor(SimuTime/deltaTau);
Spike = rand(NumTrials , NumBins) < r*deltaTau;
t = 0:deltaTau:SimuTime-deltaTau;
end

function [] = RasterPlot(spike_train,t)
    hold all;
    for trialCount = 1:size(spike_train,1)
         spikePos = t(spike_train(trialCount, :));
         for spikeCount = 1:length(spikePos)
             plot([spikePos(spikeCount) spikePos(spikeCount)], ...
             [trialCount-0.4 trialCount+0.4], 'k');
         end
    end
    ylim([0 size(spike_train, 1)+1]);
end

function [ISI,SpikesTime]=myISI(Spike)
NumTrials=size(Spike,1);
SpikeCount=sum(Spike,2);
SpikesTime=zeros(NumTrials,max(SpikeCount));
ISI=zeros(NumTrials,max(SpikeCount));
for i=1:NumTrials
    SpikesTime(i,[1:length(find(Spike(i,:)==1))])=find(Spike(i,:)==1);
    ISI(i,[1:length(find(Spike(i,:)==1))-1])=diff(SpikesTime(i,[1:length(find(Spike(i,:)==1))]),1);
end
ISI=reshape(ISI,[1,NumTrials*max(SpikeCount)]);
ISI=ISI(ISI>0)/1000;
end

function [RenewalSpike]=myRenewal(Spike,k)
NumTrials=size(Spike,1)
RenewalSpike=zeros(NumTrials,1000);
for i=1:NumTrials
       spikestime=find(Spike(i,:));
       Vector=[1:k:length(find(Spike(i,:))==1)];
       VectorTimes=spikestime(Vector);
       RenewalSpike(i,VectorTimes)=1;
end
end







