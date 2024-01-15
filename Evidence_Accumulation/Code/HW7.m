
%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_7




%% Part1.Q2-a
close all;clc;clear;
nTrial=10;
Time=1;B=1;sigma=1;dt=0.1;
for i=1:nTrial
[x(i,:),choice(i)]=simple_model(B,sigma,dt,Time);
end
%% Part1.Q2-b
close all;clc;clear;
rng shuffle
nTrial=20;
sigma=1;
dt=0.1;
Time=10;
Bvec=[1,10,0,0.1,-1];
Color={'r','b','g','k','c'};
figure;
for j=1:length(Bvec)
    B=Bvec(j);
    for i=1:nTrial
    [x,~]=simple_model(B,sigma,dt,Time);
    end
    plot(0:dt:Time-dt,x,'linewidth',1.5,'color',Color{1,j})
    hold on
end
xlabel('Time','interpreter','latex')
ylabel('x(t)','interpreter','latex')
legend('B=1','B=10','B=0','B=0.1','B=-1','location','northeastoutside')
print -depsc part2b.eps

%% Part1.Q3
close all;clc;clear;
TimeVec=[0.5:0.5:10];
Bvec=[0.1,0.2];sigma=1;dt=0.1;
nTrial=10000;
for b=1:length(Bvec)
    B=Bvec(b);
for j=1:length(TimeVec)
    Time=TimeVec(j);
    for i=1:nTrial
    [~,choice(i)]=simple_model(B,sigma,dt,Time);
    end
    ErrorRate(b,j)=(length(find(choice==0))/length(choice))*100;
end
end
scatter(0.5:0.5:10,ErrorRate(1,:),'b','filled')
hold on
scatter(0.5:0.5:10,ErrorRate(2,:),'k','filled')
hold on
plot(0.5:0.1:10,(normcdf(0,Bvec(1).*[0.5:0.1:10],sigma.*sqrt([0.5:0.1:10])))*100,'r','linewidth',2.5)
hold on
plot(0.5:0.1:10,(normcdf(0,Bvec(2).*[0.5:0.1:10],sigma.*sqrt([0.5:0.1:10])))*100,'g','linewidth',2.5)
xlim([0.5 10])
ylim([0,100])
xlabel('Time Interval(s)','interpreter','latex')
ylabel('Error Rate (\%)','interpreter','latex')
ylim([20 60])
legend('Simulation,B=0.1','Simulation,B=0.2','Theory,B=0.1','Theory,B=0.2','interpreter','latex')
print -depsc part3.eps

figure;
Error1=(normcdf(0,Bvec(1).*[0.5:0.1:10],sigma.*sqrt([0.5:0.1:10])))*100;
Error2=(normcdf(0,Bvec(2).*[0.5:0.1:10],sigma.*sqrt([0.5:0.1:10])))*100;
plot(0.5:0.1:10,(Error1-Error2)./Error1,'k','linewidth',2)
xlabel('Time Interval(s)','interpreter','latex')
ylabel('ratio','interpreter','latex')
text(4.85,0.18,'$\frac{ER_{B=0.1}-ER_{B=0.2}}{ER_{B=0.1}}$','Color','k','FontSize',15,'interpreter','latex')
print -depsc part3_a.eps

%% Part1.Q4
close all;clc;clear;
Time=10;B=0.1;sigma=1;dt=0.1;
nTrial=1000;
for i=1:nTrial
    [x(:,i),~]=simple_model(B,sigma,dt,Time);
end
meanSim=mean(x,2);
stdSim=std(x,[],2);
for i=1:nTrial
    plot(0:dt:Time-dt,x(:,i),'k','linewidth',0.1)
    hold on
end
plot(0:dt:Time-dt,meanSim,'b')
hold on
plot(0:dt:Time-dt,meanSim+stdSim,'b')
hold on
plot(0:dt:Time-dt,meanSim-stdSim,'b')
hold on 
plot(0:dt:Time-dt,B*[0:dt:Time-dt],'r')
hold on
plot(0:dt:Time-dt,B*[0:dt:Time-dt]+sigma*sqrt([0:dt:Time-dt]),'r')
hold on
plot(0:dt:Time-dt,B*[0:dt:Time-dt]-sigma*sqrt([0:dt:Time-dt]),'r')
xlabel('Time(s)','interpreter','latex')
ylabel('X(t)','interpreter','latex')
print -depsc part4.eps

%% Part1.Q5,6,7
close all;clc;clear;
Time=10;B=0.1;sigma=1;dt=0.1;
nTrial=10000;
x0=0;
ThetaP=1;ThetaN=-1;
for i=1:nTrial
[RT(i),choice(i)] = two_choice_trial(ThetaP,ThetaN,sigma,x0,B,dt)
end
save('RT.mat','RT');
save('choice.mat','choice');
x1=min(RT):0.1:max(RT);
h1 = hist(RT , length(x1));
x2=min(RT(find(choice==1))):0.1:max(RT(find(choice==1)));
h2=hist(RT(find(choice==1)) , length(x2));
x3=min(RT(find(choice==-1))):0.1:max(RT(find(choice==-1)));
h3=hist(RT(find(choice==-1)) , length(x3));
figure;
bar(x1,h1,'b');
hold on
bar(x2,h2,'g');
hold on
bar(x3,h3,'r');
xlabel('Reaction Time(s)','interpreter','latex')
ylabel('Count','interpreter','latex')
legend('Total','Correct','Incorrect')
print -depsc part7.eps

%% Part1.Q8,9
clear all; close all; clc;

sigma1 = 1;
sigma2 = 1;
ThU1=0.75;
ThU2=0.75;
ThL1=-0.75;
ThL2=-0.75;
B1 = 0.1;
B2 = 0.1;
dt = 0.01;
Time=3;
nTrial=10000;

for j=1:nTrial
[RT(j),choice(j),W(j),x1(j).val,x2(j).val]=race_trial(sigma1,sigma2,B1,B2,dt,Time,ThU1,ThU2,ThL1,ThL2);
end
figure;
set(gcf,'units','points','position',[0,0,1600,650])
Trials=randi(2,[1 nTrial]);
for k=1:length(Trials)
plot(0:dt:length(x1(Trials(k)).val)*dt-dt,x1(Trials(k)).val,'b');
hold on
plot(0:dt:length(x2(Trials(k)).val)*dt-dt,x2(Trials(k)).val,'r');
end
hold on
plot(0:dt:(max(RT(Trials))+4),ThU1*ones(1,length([0:dt:(max(RT(Trials))+4)])),'k','linewidth',2);
hold on
plot(0:dt:(max(RT(Trials))+4),ThL1*ones(1,length([0:dt:(max(RT(Trials))+4)])),'k','linewidth',2);
text((max(RT(Trials))+0.1)/2,3.15,'Choice 1','Color','k','FontSize',12)
text((max(RT(Trials))+0.1)/2,-3.15,'Choice 2','Color','k','FontSize',12)
ylim([-0.6 0.6])
xlim([0 max(RT(Trials))+0.25])
xlabel('Time','interpreter','latex')
ylabel('x(t)','interpreter','latex');
print -depsc part8a.eps

x1=min(RT):0.01:max(RT);
h1 = hist(RT , length(x1));
x2=min(RT(find(choice==1))):0.01:max(RT(find(choice==1)));
h2=hist(RT(find(choice==1)) , length(x2));
x3=min(RT(find(choice==2))):0.01:max(RT(find(choice==2)));
h3=hist(RT(find(choice==2)) , length(x3));
figure;
bar(x1,h1,'b');
hold on
bar(x2,h2,'g');
hold on
bar(x3,h3,'r');
xlabel('Reaction Time','interpreter','latex')
ylabel('Count','interpreter','latex')
legend('Total','Choose1','Choose2')
print -depsc part8b.eps

clear x2 x3 h2 h3
figure;
x2=min(RT(intersect(find(choice==1),find(W==1)))):0.01:max(RT(intersect(find(choice==1),find(W==1))));
h2=hist(RT(intersect(find(choice==1),find(W==1))) , length(x2));
x3=min(RT(intersect(find(choice==1),find(W==2)))):0.01:max(RT(intersect(find(choice==1),find(W==2))));
h3=hist(RT(intersect(find(choice==1),find(W==2))) , length(x3));
x4=min(RT(intersect(find(choice==2),find(W==1)))):0.01:max(RT(intersect(find(choice==2),find(W==1))));
h4=hist(RT(intersect(find(choice==2),find(W==1))) , length(x4));
x5=min(RT(intersect(find(choice==2),find(W==2)))):0.01:max(RT(intersect(find(choice==2),find(W==2))));
h5=hist(RT(intersect(find(choice==2),find(W==2))) , length(x5));
bar(x1,h1,'b');
hold on
bar(x2,h2,'g');
hold on
bar(x3,h3,'r');
hold on
bar(x4,h4,'y');
hold on
bar(x5,h5,'c');
xlabel('Reaction Time','interpreter','latex')
ylabel('Count','interpreter','latex')
legend('Total','Choose1,W=Racer1','Choose1,W=Racer2','Choose2,W=Racer1','Choose2,W=Racer2')
print -depsc part8c.eps

% Part 2
%% Part2.Q1
clear all; close all; clc;
rng shuffle

MT_p_values = [0.05;0.03];
LIP_weights = [0.08;-0.1];
LIP_threshold = 40;
[LIP_event_times,MT1_spikes,MT2_spikes,LIP_spikes]=LIP_activity(MT_p_values,LIP_weights,LIP_threshold)

Data=[LIP_spikes;MT1_spikes;MT2_spikes];
Data=logical(Data);
RasterPlot(Data,0:0.001:size(Data,2)*0.001-0.001);
set(gca,'ytick',[1 2 3])
print -depsc part2Q1.eps

%% Part2.Q2
clear all; close all; clc;
rng shuffle
MT_p_values_Mat = [0.05,0.03;0.03,0.05];
LIP1_weights = [-0.08;0.08];
LIP2_weights = [0.08;-0.08];
LIP1_threshold = 40;
LIP2_threshold = 40;
[rt,LIP1_event_times,LIP2_event_times,MT1_spikes,MT2_spikes,LIP1_spikes,LIP2_spikes,t_min_stimuli,t_max_stimuli]=LIP_activity2(MT_p_values_Mat,LIP1_weights,LIP2_weights,LIP1_threshold,LIP2_threshold)

Data=[LIP1_spikes;LIP2_spikes;MT1_spikes;MT2_spikes];
Data=logical(Data);
RasterPlot(Data,0:0.001:size(Data,2)*0.001-0.001);
set(gca,'ytick',[1 2 3 4])
patch([t_min_stimuli t_max_stimuli t_max_stimuli t_min_stimuli], [0 0 5 5], 'r','facealpha',0.3,'LineStyle','none')
print -depsc part2Q2.eps


function [x,choice]=simple_model(B,sigma,dt,Time)
x(1)=0;
for n=1:(Time/dt)-1
   dw=normrnd(0,sqrt(dt));
   x(n+1)=x(n)+B*dt+sigma*dw;
end
if x(Time/dt)>0
    choice=1;
else
    choice=0;
end
end


function [x,choice]=simple_model2(B,sigma,SP,Time)
p = normcdf(0,B*Time+SP,sigma*sqrt(Time));
if rand()>p
    choice=1;
else if rand()<p
        choice=-1;
    end
end
end

function [RT,choice] = two_choice_trial(ThetaP,ThetaN,sigma,x0,B,dt)
    t = 0;
    i=1;
    x(1)=0;
    while (x(i) <= ThetaP && x(i) >= ThetaN)
        x(i+1) = x(i) + B * dt + sigma * normrnd(0,sqrt(dt));
        t = t+dt;
        i=i+1;
    end
    RT = t;
    if (x(i) >= ThetaP)
        choice = 1;
    else
        choice = -1;
    end
end

function [RT,choice,W,x1,x2]=race_trial(sigma1,sigma2,B1,B2,dt,Time,ThU1,ThU2,ThL1,ThL2)
t=0;
x1(1)=0;
x2(1)=0;
for i=1:Time/dt-1
    x1(i+1)=x1(i)+B1*dt+sigma1*normrnd(0,sqrt(dt));
    x2(i+1)=x2(i)+B2*dt+sigma2*normrnd(0,sqrt(dt));
    t=t+dt;
    if x1(i)>ThU1
        choice=1;
        W=1;
        break;
    end
    if x2(i)>ThU2
        choice=1;
        W=2;
        break;
    end
    if x1(i)<ThL1
        choice=2;
        W=1;
        break;
    end
    if x2(i)<ThL2
        choice=2;
        W=2;
        break;
    end
    if t>=Time-dt
        if (abs(x1(i)-ThU1)<abs(x2(i)-ThU2))
            choice = 1;
            W=1;
        else
            choice = 1;
            W=2;
        end
        if ((x1(i)-ThL1)>(x2(i)-ThL2))
            choice = 2;
            W=2;
        else
            choice = 2;
            W=1;
        end
        break;
    end
end
RT=t;
end

function [LIP_event_times,MT1_spikes,MT2_spikes,LIP_spikes]=LIP_activity(MT_p_values,LIP_weights,LIP_threshold)
dt=0.001;
rate=0;
N=[0;0];
M=100;
t=0;
MT1_spikes=[];
MT2_spikes=[];
LIP_spikes=[];
LIP_event_times=[];
while (rate<LIP_threshold)
   dN = rand(2,1) < MT_p_values;
   N = N + dN;
   p_lip = sum(N.*LIP_weights);
   MT1_spikes = [MT1_spikes,dN(1)];
   MT2_spikes = [MT2_spikes,dN(2)];
   LIP_event = rand()<p_lip;
   LIP_spikes=[LIP_spikes,LIP_event];
   if (LIP_event == 1)
       LIP_event_times = [LIP_event_times,t];
   end 
   if (length(LIP_event_times)>=M)
       rate=M/(t-LIP_event_times(end-M+1));
   end
   t=t+dt;
end
end

function [rt,LIP1_event_times,LIP2_event_times,MT1_spikes,MT2_spikes,LIP1_spikes,LIP2_spikes,t_min_stimuli,t_max_stimuli]=LIP_activity2(MT_p_values_Mat,LIP1_weights,LIP2_weights,LIP1_threshold,LIP2_threshold)
dt=0.001;
rate1=0;
rate2=0;
N=[0;0];
M=100;
t=0;
MT1_spikes=[];
MT2_spikes=[];
LIP1_spikes=[];
LIP1_event_times=[];
LIP2_spikes=[];
LIP2_event_times=[];
t_min_stimuli=rand*0.1+0.1;
t_max_stimuli=t_min_stimuli+0.2;

while (rate1<LIP1_threshold && rate2<LIP2_threshold)
    if (t>t_min_stimuli && t<t_max_stimuli)
      MT_p_values=MT_p_values_Mat(2,:);
    else 
      MT_p_values=MT_p_values_Mat(1,:); 
    end
   dN = rand(2,1) < MT_p_values';
   N = N + dN;
   p_lip1 = sum(N.*LIP1_weights);
   p_lip2 = sum(N.*LIP2_weights);
   MT1_spikes = [MT1_spikes,dN(1)];
   MT2_spikes = [MT2_spikes,dN(2)];
   LIP1_event = rand()<p_lip1;
   LIP2_event = rand()<p_lip2;
   LIP1_spikes=[LIP1_spikes,LIP1_event];
   LIP2_spikes=[LIP2_spikes,LIP2_event];
   if (LIP1_event == 1)
       LIP1_event_times = [LIP1_event_times,t];
   end 
   if (LIP2_event == 1)
       LIP2_event_times = [LIP2_event_times,t];
   end 
   if (length(LIP1_event_times)>=M)
       rate1=M/(t-LIP1_event_times(end-M+1));
   end
   if (length(LIP2_event_times)>=M)
       rate2=M/(t-LIP2_event_times(end-M+1));
   end
   t=t+dt;
end
rt=t;
end

function [] = RasterPlot(spike_train,t)
    hold all;
    for trialCount = 1:size(spike_train,1)
         spikePos = t(spike_train(trialCount, :));
         for spikeCount = 1:length(spikePos)
             plot([spikePos(spikeCount) spikePos(spikeCount)], ...
             [trialCount-0.1 trialCount+0.1], 'k');
         end
    end
    ylim([0 size(spike_train, 1)+1]);
    xlabel('Time(s)','interpreter','latex')
end

