%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_1

%% Part2 
%% section a


clear;
clc;
close all;

SimulationTime=0.1;
deltaT=0.0001;
tauM=0.013;
t=0:deltaT:SimulationTime;
RI=0.02*ones(1,length(t));
RefPeriod=0;
[v,nSpike]=myLIF(RI,SimulationTime,tauM,t,RefPeriod,deltaT)

figure;
plot(t,v,'b');
title('LIF Model with constant input current of 20 mv','interpreter', 'latex');
xlabel('Time (s)','interpreter', 'latex');
ylabel('Voltage (mv)','interpreter', 'latex');
print -depsc fig9.eps

%% section c
r=300; %% frequency 
SimulationTime=1; 
tauM=0.013;
deltaT=0.0001;
t = 0:deltaT:SimulationTime-deltaT;
NumTrials=1;
Spike=SpikePoissonGen(r,SimulationTime,deltaT,NumTrials);
t_peak=0.0015;
tKernel=0:deltaT:0.05-deltaT;
kernel=tKernel.*exp(-tKernel./t_peak);
I0=0.1;
I=conv(Spike,kernel,'same');
I=I/max(I);
I=I0*I;

figure;
set(gcf,'units','points','position',[0,0,1000,600])
subplot(3,1,1)
hold all;
spikePos = t(Spike(1, :));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
           [0 1], 'k','color','b');
         end
xlabel('Time (s)','interpreter', 'latex');
ylabel('Spike','interpreter', 'latex');
title(['Input Spike Train'],'fontsize',8,'interpreter', 'latex')
ylim([0 2])
subplot(3,1,2)
plot(t,I,'color','b')
xlabel('Time (s)','interpreter', 'latex');
ylabel('I(t) [A]','interpreter', 'latex');
title(['Realistic I(t)'],'fontsize',8,'interpreter', 'latex')
ylim([0 0.12])
ref=0.001;
[v,nSpike]=myLIF(I,SimulationTime,tauM,t,ref,deltaT)
subplot(3,1,3)
plot(t,v,'b');
xlabel('Time (s)','interpreter', 'latex');
ylabel('V(v)','interpreter', 'latex');
title(['Realistic I(t)'],'fontsize',8,'interpreter', 'latex')
ylim([0 0.031])
hold on
plot(t, 0.015*ones(1,length(t)),'color','r','LineStyle','--');
print -depsc fig10.eps

% The effect of width and magnitude of EPSCs on CV
I0=0.05:0.001:0.3;
I=conv(Spike,kernel,'same');
I=I/max(I);
for i=1:length(I0)
    I_i=I0(i)*I;
   [v,nSpike]=myLIF(I_i,SimulationTime,tauM,t,ref)
    ISI=diff(find(v==0.030))*deltaT/1000;
    cv(i)=std(ISI)/mean(ISI);
end
figure;
scatter(I0,cv)
xlabel('$I_0$','interpreter', 'latex');
ylabel('CV','interpreter', 'latex');
title(['Relationship between $I_0$ and CV'],'fontsize',8,'interpreter', 'latex')
print -depsc fig11.eps
 
t_peak=0.001:0.001:0.2; 
for i=1:length(t_peak)
    kernel=tKernel.*exp(-tKernel./t_peak(i));
    I=conv(Spike,kernel,'same');
    I=I/max(I);
    I0=0.1;
    I=I0*I;
    [v,nSpike]=myLIF(I,SimulationTime,tauM,t,ref)
    ISI=diff(find(v==0.030))*deltaT/1000;
    cv(i)=std(ISI)/mean(ISI);
end
scatter(t_peak,cv)
xlabel('$t_{peak}$','interpreter', 'latex');
ylabel('CV','interpreter', 'latex');
title(['Relationship between $t_{peak}$ and CV'],'fontsize',8,'interpreter', 'latex')
print -depsc fig12.eps

%% section d
r=50;
SimulationTime=1; 
tauM=0.013;
deltaT=0.0001;
t = 0:deltaT:SimulationTime-deltaT;
t_peak=0.0015;
tKernel=0:deltaT:0.05-deltaT;
NumTrials=1;
ref=0.001;
I0=0.2;
NumNeuron=100;
p=[0:0.5:35];
for i=1:length(p)
    all=NeuronSelect(NumNeuron,p(i));
for j=1:length(all)
    if all(j)==1
            Spike=SpikePoissonGen(r,SimulationTime,deltaT,NumTrials);
    else if all(j)==0
            Spike=SpikePoissonGen(r,SimulationTime,deltaT,NumTrials);
            Spike=-Spike;
        end
    end
    kernel=tKernel.*exp(-tKernel./t_peak);
    I(j,:)=conv(Spike,kernel,'same');
    I(j,:)=I0*(I(j,:)/max(abs(I)));
    end
    I=sum(I,1);
    [v, nSpike]=myLIF(I,SimulationTime,tauM,t,ref,deltaT)
    z=diff(find(v==0.03))/10000;
    cv(i)=std(z)/mean(z);
end
scatter(p,cv)
xlabel('Inhibition percentage of synaptic inputs','interpreter', 'latex');
ylabel('CV','interpreter', 'latex');
print -depsc fig13.eps

%% section e
D=0.02:0.001:0.1;
NoverM=0.1:0.025:0.8;
CV=zeros(length(D),length(NoverM))
M=250; 
r=50 % freq
SimuTime=1;
deltaTau=0.0001;
t=0:deltaTau:1-deltaTau;

Spike=SpikePoissonGen(r,SimuTime,deltaTau,M);
for j=1:length(NoverM)
    N=NoverM(j)*M;
    Spike_new=zeros(1,length(t));
for i=1:length(D)
    [nSpikeExc,nSpikeInh,Nwidth]=nSpikeD_Interval(Spike,D(i),deltaTau,t)
    for k=1:size(nSpikeExc,2)
    if (nSpikeExc(k)>=N)
        idx=k*Nwidth;
        if idx>=length(t)
            idx=length(t)
        end
         Spike_new(idx)=1
    end
    end
    ISI=diff(find(Spike_new==1))*deltaTau;
    CV(i,j)=std(ISI)/mean(ISI);
end
end
figure;
set(gcf,'units','points','position',[0,0,500,400])
for i=1:length(D)
    scatter(D,CV(:,i))
    hold on
    xlabel('D(s)','interpreter', 'latex');
    ylabel('CV','interpreter', 'latex');
    title(['Relationship between CV and D'],'fontsize',8,'interpreter', 'latex')
    ylim([0 2])
end
print -depsc fig13.eps

figure;
set(gcf,'units','points','position',[0,0,500,400])
for i=1:length(NoverM)
   scatter(NoverM,CV(i,:))
    hold on
    xlabel('$\frac{N}{M}$','interpreter', 'latex');
    ylabel('CV','interpreter', 'latex');
    title(['Relationship between CV and $\frac{N}{M}$'],'fontsize',8,'interpreter', 'latex')
    ylim([0 2])
end
 print -depsc fig14.eps

%% section f

   allNeurons=250;
   excitatoryN=200;
   inhibN=50;
   dt=0.0001;
   SimuTime=1;
   r=50; %freq
   t=0:dt:SimuTime-dt;
   inhibNeurons=randi([1 allNeurons],1,inhibN);
   Spike=SpikePoissonGen(r,SimuTime,dt,allNeurons);
   Spike_new=zeros(allNeurons,length(t));
   for i=1:length(inhibNeurons)
       Spike_new(inhibNeurons(i),:)=Spike(inhibNeurons(i),:).*(-2);
   end
   Spike=Spike+Spike_new;
   Thre=5:5:125;
   D=0.001:0.001:0.1;
   
   for i=1:length(Thre)
       spike_new=zeros(1,length(t));
       for j=1:length(D)
          [nSpikeExc,nSpikeInh,Nwidth]=nSpikeD_Interval(Spike,D(j),dt,t) 
             for k=1:size(nSpikeExc,2)
    if (nSpikeExc(k)-nSpikeInh(k)>=Thre(i))
        idx=k*Nwidth;
        if idx>=length(t)
            idx=length(t)
        end
         spike_new(idx)=1
    end
    end
    ISI=diff(find(spike_new==1))*dt;
    CV(i,j)=std(ISI)/mean(ISI);
       end
   end
  figure;
  set(gcf,'units','points','position',[0,0,450,350])
for i=1:length(D)
    scatter(D,CV(i,:))
    hold on
    xlabel('D(s)','interpreter', 'latex');
    ylabel('CV','interpreter', 'latex');
    title(['Relationship between CV and D'],'fontsize',8,'interpreter', 'latex')
    ylim([0 1.5])
end
print -depsc fig15.eps

figure;
set(gcf,'units','points','position',[0,0,450,350])
for i=1:length(Thre)
    scatter(Thre,CV(:,i))
    hold on
    xlabel('Threshold','interpreter', 'latex');
    ylabel('CV','interpreter', 'latex');
    title(['Relationship between CV and Threshold'],'fontsize',8,'interpreter', 'latex')
    ylim([0 2])
end
   print -depsc fig16.eps


%% functions

function [v,nSpike]=myLIF(I,SimulationTime,tauM,t,RefPeriod,deltaT)
Vr=0;
Vpeak=0.03;
Vth=0.015;
v=zeros(1,length(t));
nSpike=0;
last_spike_time=-inf;
for i=2:length(t)
    if v(i-1)>Vth && (t(i)-last_spike_time)>RefPeriod
        v(i-1)=Vpeak;
        v(i)=Vr
        nSpike=nSpike+1;
        last_spike_time=t(i);
    else
        dv=(1/tauM)*(-v(i-1)+I(i))*deltaT;
        v(i)=v(i-1)+dv;
    end
end
end


function Spike=SpikePoissonGen(r,SimuTime,deltaTau,NumTrials)
NumBins=floor(SimuTime/deltaTau);
Spike = rand(NumTrials , NumBins) < r*deltaTau;
end

function [all]=NeuronSelect(NumNeuron,p)
NumInhib=floor((p/100)*NumNeuron);
inhibitory=zeros(1,NumInhib);
NumExcit=floor(((100-p)/100)*NumNeuron);
excitatory=ones(1,NumExcit);
all=[inhibitory excitatory];
all=all(randperm(length(all)));
end

function [nSpikeExc,nSpikeInh,Nwidth]=nSpikeD_Interval(Spike,D,dt,t)
Nt=length(t);
Nwidth = round(D(1)/dt);       
Nwindows = round(Nt/Nwidth);         
nSpikeExc = zeros(1,Nwindows);  
nSpikeInh = zeros(1,Nwindows);  
for i = 1:Nwindows
    idx=i*Nwidth;
    if idx>=Nt
        idx=Nt;
    end
    nSpikeExc(i) = length(find(sum(Spike(:,(i-1)*Nwidth+1:idx),2)>0));
    nSpikeInh(i) = length(find(sum(Spike(:,(i-1)*Nwidth+1:idx),2)<0));
end
end