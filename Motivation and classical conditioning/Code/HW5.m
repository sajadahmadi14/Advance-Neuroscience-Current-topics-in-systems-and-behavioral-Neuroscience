%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_5


%% Part1
close all;clc;clear;
epsilon=0.03;
n=400;
%extinction
r=[ones(1,n/2) zeros(1,n/2)];
u=ones(1,n);
stimuliNum=1;
w=zeros(stimuliNum,n);
for i=1:n-1
    w(:,i+1)=w(:,i)+epsilon*(r(i)-w(:,i)'*u(:,i))*u(:,i);
end
figure;
set(gcf,'units','points','position',[0,0,1200,600])
subplot(2,3,1)
scatter([1:n],w,'b','filled')
title('Extinction','interpreter','latex');
xlabel('Trials','interpreter','latex');
ylabel('w','interpreter','latex');
%partial
alpha=[0.8,0.4];
u=ones(1,n);
stimuliNum=1;
w0=zeros(stimuliNum,1);
for j=1:length(alpha)
r=rand(1,n)<alpha(j);
w=zeros(stimuliNum,n);
for i=1:n-1
    w(:,i+1)=w(:,i)+epsilon*(r(i)-w(:,i)'*u(:,i))*u(:,i);
end
subplot(2,3,2)
color=['b','r']
scatter([1:n],w,color(j),'filled')
hold on
legend(['\alpha=',num2str(alpha(1))],['\alpha=',num2str(alpha(2))],'location','northeastoutside','interpreter','latex');
title('Partial','interpreter','latex');
xlabel('Trials','interpreter','latex');
ylabel('w','interpreter','latex');
end
%blocking
clear w
stimuliNum=2;
w0=zeros(stimuliNum,1);
r=[1,ones(1,n)];
if stimuliNum==2;
    u1=ones(1,n);
    u2=[zeros(1,n/2) ones(1,n/2)];
    u=[u1;u2];
end
w=zeros(stimuliNum,n);
for i=1:n-1
    w(:,i+1)=w(:,i)+epsilon*(r(i)-w(:,i)'*u(:,i))*u(:,i);
end
subplot(2,3,3)
color=['b','r']
scatter([1:n],w(1,:),color(1),'filled')
hold on
scatter([1:n],w(2,:),color(2),'filled')
legend('w_1','w_2','location','northeastoutside','interpreter','latex');
title('blocking','interpreter','latex');
xlabel('Trials','interpreter','latex');
ylabel('w','interpreter','latex');
%inhibitory
clear w
stimuliNum=2;
if stimuliNum==2;
    u1=ones(1,n);
    u2=rand(1,n)<0.3;
    u=[u1;u2];
end
r=~u2;
w=zeros(stimuliNum,n);
for i=1:n-1
    w(:,i+1)=w(:,i)+epsilon*(r(i)-w(:,i)'*u(:,i))*u(:,i);
end
subplot(2,3,4)
color=['b','r']
scatter([1:n],w(1,:),color(1),'filled')
hold on
scatter([1:n],w(2,:),color(2),'filled')
legend('w_1','w_2','location','northeastoutside','interpreter','latex');
title('Inhibitory','interpreter','latex');
xlabel('Trials','interpreter','latex');
ylabel('w','interpreter','latex');
%overshadow
clear w
stimuliNum=2;
r=[1,ones(1,n)];
if stimuliNum==2;
    u1=ones(1,n);
    u2=rand(1,n)<1;
    u=[u1;u2];
end
r=ones(1,n);
w=zeros(stimuliNum,n);
for i=1:n-1
    w(:,i+1)=w(:,i)+epsilon*(r(i)-w(:,i)'*u(:,i))*u(:,i);
end
subplot(2,3,5)
color=['b','r']
scatter([1:n],w(1,:),color(1),'filled')
hold on
scatter([1:n],w(2,:),color(2),'filled')
legend('w_1','w_2','location','northeastoutside','interpreter','latex');
title('Overshadow','interpreter','latex');
xlabel('Trials','interpreter','latex');
ylabel('w','interpreter','latex');
print -depsc part1.eps

% Part2
%% blocking and unblocking , Drift (Figure 1B,C,D)
close all;clc;clear;
n=20;
r=ones(1,n)
xMat=[ones(1,n);zeros(1,n/2) ones(1,n/2)];
stimuliNum=2;
tau=0.5;
processNoise=0.003;
sigma0=0.7;
[w,sigma1,sigma2]=KalmanFilter(xMat,n,r,stimuliNum,tau,processNoise,sigma0)
figure;
set(gcf,'units','points','position',[0,0,1200,400])
subplot(2,2,1)
plot(1:n,w(1,:),'k','linewidth',2)
hold on 
plot(n/2:n,w(2,n/2:n),'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('Trials','interpreter','latex')
ylabel('w(t)','interpreter','latex')
title('Blocking ($\tau=0.5$ , $\sigma_0=0.7$ , ProcessNoise=0.003)','interpreter','latex')
set(gca,'ytick',[0:0.1:1])
subplot(2,2,3)
plot(1:n,sigma1,'k','linewidth',2)
hold on 
plot(n/2:n-1,sigma2,'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('Trials','interpreter','latex')
ylabel('$\sigma^2(t)$','interpreter','latex')
title('Blocking','interpreter','latex')
set(gca,'ytick',[0:0.1:1])
r=[ones(1,n/2),2*ones(1,n/2)];
[w,sigma1,sigma2]=KalmanFilter(xMat,n,r,stimuliNum,tau,processNoise,sigma0)
subplot(2,2,2)
plot(1:n,w(1,:),'k','linewidth',2)
hold on 
plot(n/2:n,w(2,n/2:n),'--k','linewidth',2)
xlim([1 n]);
ylim([0 1.2]);
xlabel('Trials','interpreter','latex')
ylabel('w(t)','interpreter','latex')
title('Unblocking ($\tau=0.5$ , $\sigma_0=0.7$ , ProcessNoise=0.003)','interpreter','latex')
set(gca,'ytick',[0:0.1:1])
subplot(2,2,4)
plot(1:n,sigma1,'k','linewidth',2)
hold on 
plot(n/2:n-1,sigma2,'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('Trials','interpreter','latex')
ylabel('$\sigma^2(t)$','interpreter','latex')
title('Unblocking','interpreter','latex')
set(gca,'ytick',[0:0.1:1])
print -depsc blocking_unblocking.eps

clear w
rng shuffle 
w=zeros(2,n);
w(:,1)=1;
v=zeros(2,n);
v(1,:)=normrnd(0,0.1,1,n);
v(2,:)=normrnd(0,0.1,1,n);
for i=1:n-1
    w(1,i+1)=w(1,i)+v(1,i);
    w(2,i+1)=w(2,i)+v(2,i);
end
figure;
plot(1:n,w(1,:),'k');
hold on 
plot(1:n,w(2,:),'--k');
xlabel('Trials','interpreter','latex')
ylabel('w','interpreter','latex')
xlim([1 20]);
ylim([0 2]);
legend('w_1(t)','w_2(t)','interpreter','latex')
print -depsc drift.eps

%% Backwards blocking (Figure2)
close all;clc;clear;
n=20;
r=ones(1,n)
xMat=[ones(1,n);ones(1,n/2),zeros(1,n/2)];
stimuliNum=2;
tau=0.9;
processNoise=0.03;
sigma0=0.7;
[w,sigma1,sigma2,sigmaMat]=KalmanFilter(xMat,n,r,stimuliNum,tau,processNoise,sigma0)
sigmaMatP = {sigmaMat{1}, sigmaMat{5},sigmaMat{9}, sigmaMat{13}, sigmaMat{16},sigmaMat{19}};
wMat = {w(:,1), w(:,5), w(:,9) , w(:,13),w(:,16), w(:,19)};
r = 2:-0.2:0;
colorCode = linspace(0,1,length(r));
for i = 1:length(colorCode)
    color{i} = [colorCode(i),colorCode(i),colorCode(i)];
end
theta = 0:0.1:2*pi;
figure
for sub = 1:6
    subplot(2,3,sub)
    tmpW = wMat{sub};
    sigmaTmp = sigmaMatP{sub};
    for i = 1:length(r)
        x = r(i)*sin(theta);
        y = r(i)*cos(theta);
        for j = 1:length(x)
            tmp = [x(j);y(j)];
            tmp = sigmaTmp * tmp;
            x(j) = tmp(1)+tmpW(1);
            y(j) = tmp(2)+tmpW(2);
        end
        fill(x,y,color{i},'LineStyle','none')
        axis square
        hold on
    end
    scatter(x(j),y(j),'*k')
    t = {'1' '5' '9' '13' '16' '19'};
    xlim([-1 2])
    ylim([-1 2])
    set(gca,'Color','k')
    xlabel('$\bar{\omega}_1$','interpreter','LaTex')
    ylabel('$\bar{\omega}_2$','interpreter','LaTex')
    title(['t= ', t{sub}],'interpreter','latex')
end
print -depsc fig2Paper.eps
%% Changing process noise (Blocking and Unblocking)
close all;clc;clear;
n=20;
r=ones(1,n)
xMat=[ones(1,n);zeros(1,n/2) ones(1,n/2)];
stimuliNum=2;
tau=0.5;
processNoiseVec=[0.001,0.005,0.01,0.03];
ProNoiVec={'0.001','0.005','0.01','0.03'};
sigma0=0.7;
figure;
set(gcf,'units','points','position',[0,0,1200,400])
for i=1:length(processNoiseVec)
[w,sigma1,sigma2]=KalmanFilter(xMat,n,r,stimuliNum,tau,processNoiseVec(i),sigma0)
subplot(2,4,i)
plot(1:n,w(1,:),'k','linewidth',2)
hold on 
plot(n/2:n,w(2,n/2:n),'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('Trials','interpreter','latex')
ylabel('w(t)','interpreter','latex')
title(['Blocking ($\tau=0.5$ , $\sigma_0=0.7$ , ProcessNoise=',ProNoiVec{1,i},')'],'interpreter','latex','fontsize',8)
%set(gca,'xtick',[1:2:18])
set(gca,'ytick',[0:0.1:1])
subplot(2,4,i+4)
plot(1:n,sigma1,'k','linewidth',2)
hold on 
plot(n/2:n-1,sigma2,'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('Trials','interpreter','latex')
ylabel('$\sigma^2(t)$','interpreter','latex')
title(['Blocking ($\tau=0.5$ , $\sigma_0=0.7$ , ProcessNoise=',ProNoiVec{1,i},')'],'interpreter','latex','fontsize',8)
%set(gca,'xtick',[0:2:19])
set(gca,'ytick',[0:0.1:1])
end
print -depsc changProcess1_block.eps

r=[ones(1,n/2),2*ones(1,n/2)];
figure;
set(gcf,'units','points','position',[0,0,1200,400])
for i=1:length(processNoiseVec)
[w,sigma1,sigma2]=KalmanFilter(xMat,n,r,stimuliNum,tau,processNoiseVec(i),sigma0)
subplot(2,4,i)
plot(1:n,w(1,:),'k','linewidth',2)
hold on 
plot(n/2:n-1,w(2,n/2:n-1),'--k','linewidth',2)
xlim([1 n]);
ylim([0 1.1]);
xlabel('Trials','interpreter','latex')
ylabel('w(t)','interpreter','latex')
title(['Unblocking ($\tau=0.5$ , $\sigma_0=0.7$ , ProcessNoise=',ProNoiVec{1,i},')'],'interpreter','latex','fontsize',8)
%set(gca,'xtick',[0:2:18])
set(gca,'ytick',[0:0.1:1])
subplot(2,4,i+4)
plot(1:n,sigma1,'k','linewidth',2)
hold on 
plot(n/2:n-1,sigma2,'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('Trials','interpreter','latex')
ylabel('$\sigma^2(t)$','interpreter','latex')
title(['Unblocking ($\tau=0.5$ , $\sigma_0=0.7$ , ProcessNoise=',ProNoiVec{1,i},')'],'interpreter','latex','fontsize',8)
%set(gca,'xtick',[0:2:19])
set(gca,'ytick',[0:0.1:1])
end
print -depsc changProcess1_unblock.eps
%% Changing process noise (Backwards Blocking)
clear all
processNoiseVec=[0.01,0.03,0.06,0.09];
for k=1:length(processNoiseVec)
n=20;
r=ones(1,n);
xMat=[ones(1,n);zeros(1,n/2) ones(1,n/2)];
stimuliNum=2;
ProNoiVec={'0.01','0.03','0.06','0.09'};
stimuliNum=2;
tau=0.9;
sigma0=0.7;
[w,sigma1,sigma2,sigmaMat]=KalmanFilter(xMat,n,r,stimuliNum,tau,processNoiseVec(k),sigma0)
sigmaMatP = {sigmaMat{1},sigmaMat{9}, sigmaMat{19}};
wMat = {w(:,1), w(:,9) , w(:,19)};
r = 2:-0.2:0;
colorCode = linspace(0,1,length(r));
for i = 1:length(colorCode)
    color{i} = [colorCode(i),colorCode(i),colorCode(i)];
end
theta = 0:0.1:2*pi;
f=figure
set(gcf,'units','points','position',[0,0,1200,400])
for sub = 1:3
    subplot(1,3,sub)
    tmpW = wMat{sub};
    sigmaTmp = sigmaMatP{sub};
    for i = 1:length(r)
        x = r(i)*sin(theta);
        y = r(i)*cos(theta);
        for j = 1:length(x)
            tmp = [x(j);y(j)];
            tmp = sigmaTmp * tmp;
            x(j) = tmp(1)+tmpW(1);
            y(j) = tmp(2)+tmpW(2);
        end
        fill(x,y,color{i},'LineStyle','none')
        axis square
        hold on
    end
    scatter(x(j),y(j),'*k')
    t = {'1' '9' '19'};
    xlim([-1 2])
    ylim([-1 2])
    set(gca,'Color','k')
    xlabel('$\bar{\omega}_1$','interpreter','LaTex')
    ylabel('$\bar{\omega}_2$','interpreter','LaTex')
    title(['t= ', t{sub},'  , ProcessNoise=',ProNoiVec{k}],'interpreter','latex')
end
    ProNoiVec={'1','3','6','9'};
    fig_name = strcat('paper1ProcessBlock_',ProNoiVec{k});
    print(f,fig_name,'-depsc'); 
end
%% Changing measurement noise (Blocking and Unblocking)
close all;clc;clear;
n=20;
r=ones(1,n)
xMat=[ones(1,n);zeros(1,n/2) ones(1,n/2)];
stimuliNum=2;
tauVec=[0.2,0.5,0.8,1.1];
TauVec={'0.2','0.5','0.8','1.1'};
processNoise=0.003;
sigma0=0.7;
figure;
set(gcf,'units','points','position',[0,0,1200,400])
for i=1:length(tauVec)
[w,sigma1,sigma2]=KalmanFilter(xMat,n,r,stimuliNum,tauVec(i),processNoise,sigma0)
subplot(2,4,i)
plot(1:n,w(1,:),'k','linewidth',2)
hold on 
plot(n/2:n,w(2,n/2:n),'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('t','interpreter','latex')
ylabel('w(t)','interpreter','latex')
title(['Blocking ($\tau=$',TauVec{1,i},' , $\sigma_0=0.7$ , ProcessNoise=0.003)'],'interpreter','latex','fontsize',8)
set(gca,'ytick',[0:0.1:1])
subplot(2,4,i+4)
plot(1:n,sigma1,'k','linewidth',2)
hold on 
plot(n/2:n-1,sigma2,'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('t','interpreter','latex')
ylabel('$\sigma^2(t)$','interpreter','latex')
title(['Blocking ($\tau=$',TauVec{1,i},' , $\sigma_0=0.7$ , ProcessNoise=0.003)'],'interpreter','latex','fontsize',8)
set(gca,'ytick',[0:0.1:1])
end
print -depsc changMesu1_block.eps

r=[ones(1,n/2),2*ones(1,n/2)];
figure;
set(gcf,'units','points','position',[0,0,1200,400])
for i=1:length(tauVec)
[w,sigma1,sigma2]=KalmanFilter(xMat,n,r,stimuliNum,tauVec(i),processNoise,sigma0)
subplot(2,4,i)
plot(1:n,w(1,:),'k','linewidth',2)
hold on 
plot(n/2:n,w(2,n/2:n),'--k','linewidth',2)
xlim([1 n]);
ylim([0 1.1]);
xlabel('t','interpreter','latex')
ylabel('w(t)','interpreter','latex')
title(['Unblocking ($\tau=$',TauVec{1,i},' , $\sigma_0=0.7$ , ProcessNoise=0.003)'],'interpreter','latex','fontsize',8)
set(gca,'ytick',[0:0.1:1])
subplot(2,4,i+4)
plot(1:n,sigma1,'k','linewidth',2)
hold on 
plot(n/2:n-1,sigma2,'--k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('t','interpreter','latex')
ylabel('$\sigma^2(t)$','interpreter','latex')
title(['Unblocking ($\tau=$',TauVec{1,i},' , $\sigma_0=0.7$ , ProcessNoise=0.003)'],'interpreter','latex','fontsize',8)
set(gca,'ytick',[0:0.1:1])
end
print -depsc changMesu1_unblock.eps

%% Changing measurement noise (Backwards Blocking)
clear all
tauVec=[0.6,0.8,1,1.2];
for k=1:length(tauVec)
n=20;
r=ones(1,n);
xMat=[ones(1,n);zeros(1,n/2) ones(1,n/2)];
stimuliNum=2;
sigma0=0.7;
processNoise=0.03;
TauVec={'0.6','0.8','1','1.2'};
[w,sigma1,sigma2,sigmaMat]=KalmanFilter(xMat,n,r,stimuliNum,tauVec(k),processNoise,sigma0)
sigmaMatP = {sigmaMat{1},sigmaMat{9}, sigmaMat{19}};
wMat = {w(:,1), w(:,9) , w(:,19)};
r = 2:-0.2:0;
colorCode = linspace(0,1,length(r));
for i = 1:length(colorCode)
    color{i} = [colorCode(i),colorCode(i),colorCode(i)];
end
theta = 0:0.1:2*pi;
f=figure
set(gcf,'units','points','position',[0,0,1200,400])
for sub = 1:3
    subplot(1,3,sub)
    tmpW = wMat{sub};
    sigmaTmp = sigmaMatP{sub};
    for i = 1:length(r)
        x = r(i)*sin(theta);
        y = r(i)*cos(theta);
        for j = 1:length(x)
            tmp = [x(j);y(j)];
            tmp = sigmaTmp * tmp;
            x(j) = tmp(1)+tmpW(1);
            y(j) = tmp(2)+tmpW(2);
        end
        fill(x,y,color{i},'LineStyle','none')
        axis square
        hold on
    end
    scatter(x(j),y(j),'*k')
    t = {'1' '9' '19'};
    xlim([-1 2])
    ylim([-1 2])
    set(gca,'Color','k')
    xlabel('$\bar{\omega}_1$','interpreter','LaTex')
    ylabel('$\bar{\omega}_2$','interpreter','LaTex')
    title(['t= ', t{sub},'  , $\tau$=',TauVec{k}],'interpreter','latex')
end
    TauVec={'6','8','10','12'};
    fig_name = strcat('paper1Tau_',TauVec{k});
    print(f,fig_name,'-depsc'); 
end
%% s1->r and s2->-r paradigm
close all;clc;clear;
n=20;
r=[ones(1,n/2),-ones(1,n/2)];
xMat=ones(1,n);
stimuliNum=1;
tau=0.7;
processNoise=0.03;
sigma0=0.7;
w=zeros(stimuliNum,n);
sigmaMat=cell(1,n);
% stimate 
sigmaMat{1}=sigma0;
sigma=[];
for i=2:n
    sigma_stimate=sigmaMat{i-1}+processNoise;
    w_stimate=w(:,i-1);
    sigmaMat{i}=sigma_stimate-sigma_stimate*xMat(:,i)*(xMat(:,i)'*sigma_stimate*xMat(:,i)+tau^2)^-1*xMat(:,i)'*sigma_stimate;
    w(:,i)=w_stimate+sigma_stimate*xMat(:,i)*(xMat(:,i)'*sigma_stimate*xMat(:,i)+tau^2)^-1*(r(i)-xMat(:,i)'*w_stimate);
    sigma=[sigma sigmaMat{1,i}(1,1)];
end
sigma=[sigma0 sigma];
figure;
set(gcf,'units','points','position',[0,0,1200,400])
subplot(1,2,1)
plot(1:n,w(1,:),'k','linewidth',2)
xlim([1 n]);
ylim([-1 1]);
xlabel('Trials','interpreter','latex')
ylabel('w(t)','interpreter','latex')
title('$\tau=0.7$ , $\sigma_0=0.7$ , ProcessNoise=0.03','interpreter','latex')
subplot(1,2,2)
plot(1:n,sigma,'k','linewidth',2)
xlim([1 n]);
ylim([0 1]);
xlabel('Trials','interpreter','latex')
ylabel('$\sigma^2(t)$','interpreter','latex')
title('$\tau=0.7$ , $\sigma_0=0.7$ , ProcessNoise=0.03','interpreter','latex')
print -depsc Q5.eps
%% Fig 3A,B,C
close all;clc;clear;
rng shuffle;
n=100;
w=zeros(1,n);
processNoise=0.1;
w(1)=normrnd(0,processNoise,1);
v=normrnd(0,processNoise,1,n);
phi=normrnd(0,1.5,1,n);
phi(46)=-2;phi(93)=4; 
c=[zeros(1,45),1,zeros(1,46),1,zeros(1,7)];
for i = 2:n
    w(i) = w(i-1) + v(i-1) + c(i-1) * phi(i-1);
end
v=normrnd(0,0.5,1,n); 
for i = 2:n
    r(i) = w(i-1) + v(i-1) + c(i-1) * phi(i-1);
end
figure;
set(gcf,'units','points','position',[0,0,800,155])
scatter(0:n-1,w,'filled','k')
hold on
scatter(0:n-1,r,'b','x')
hold on
xlim([0 n]);
ylim([-4 4]);
set(gca,'xtick',[0,46,93,100])
%set(gca,'ytick',[-4,-2,-2,0,2,3,4])
legend('w(t)','r(t)','interpreter','latex','location','northeastoutside')
print -depsc 3A.eps
figure;
set(gcf,'units','points','position',[0,0,800,155])
xMat=[ones(1,n)];
stimuliNum=1;
tau=0.7;
processNoise=0.1;
sigma0=0.7;
gamma=3.3;
[what,sigma,beta]=KalmanFilter2(xMat,n,r,stimuliNum,tau,processNoise,sigma0,gamma)
scatter(0:n-1,r,'b','x')
hold on
scatter(0:n-1,what,'r');
hold on 
scatter(0:n-1,w,'filled','k');
xlim([0 n])
ylim([-4 4])
legend('r(t)','\hat{w}(t)','w(t)','interpreter','latex','location','northeastoutside')
print -depsc 3B.eps
figure;
set(gcf,'units','points','position',[0,0,800,155])
plot(sigma,'--k','LineWidth',1.5)
hold on
plot(beta,'k','LineWidth',1.5)
hold on
plot(0:n-1,3.3*ones(1,n),'-.','color','r');
legend('ACh','NE','\gamma','Interpreter','LaTex','location','northeastoutside')
set(gca,'ytick',[0,3.3,5,10])
ylim([0 10])
print -depsc 3C.eps
%% Fig 3D
close all;clc;clear;
MSE=zeros(1000,26);
stimuliNum=1;
tau=0.6;
processNoise=0.01;
sigma0=0.6;
gamma=0:2:50;
for iter=1:2000
    n=100;
    xMat=[ones(1,n)];
w=zeros(1,n);
rng shuffle;
w(1)=normrnd(0,0.1,1);
v=normrnd(0,0.1,1,n);
phi=normrnd(0,1.5,1,n);
phi(46)=-2;phi(93)=4; 
c=[zeros(1,45),1,zeros(1,46),1,zeros(1,7)];
r = w;
for i = 2:n
    w(i) = w(i-1) + v(i-1) + c(i-1) * phi(i-1);
end
v=normrnd(0,0.5,1,n); 
r(1) = w(1) + v(1) + c(1)*phi(1);
for i = 2:n
    r(i) = w(i) + v(i) + c(i) * phi(i);
end
for j=1:length(gamma)
[what,sigma,beta]=KalmanFilter2(xMat,n,r,stimuliNum,tau,processNoise,sigma0,gamma(j))
MSE(iter,j)=sum((what-w).^2);  
end
end
errorMSE=[];
MeanMSE=[];
for i=1:26
errorMSE=[errorMSE std(MSE(:,i))];
MeanMSE=[MeanMSE median(MSE(:,i))];
end
figure;
errorbar(gamma,MeanMSE,errorMSE/10,'linewidth',1.25)
xlabel('$\gamma','interpreter','latex')
ylabel('MSE','interpreter','latex')
print -depsc 3D.eps

%% functions
function [w,sigma1,sigma2,sigmaMat]=KalmanFilter(xMat,n,r,stimuliNum,tau,processNoise,sigma0)
w=zeros(stimuliNum,n);
sigmaMat=cell(1,n);
sigmaMat{1}=[sigma0 0;0 sigma0];
sigma1=[];
sigma2=[];
for i=2:n
    sigma_stimate=sigmaMat{i-1}+[processNoise,0;0,processNoise];
    w_stimate=w(:,i-1);
    sigmaMat{i}=sigma_stimate-sigma_stimate*xMat(:,i)*(xMat(:,i)'*sigma_stimate*xMat(:,i)+tau^2)^-1*xMat(:,i)'*sigma_stimate;
    w(:,i)=w_stimate+sigma_stimate*xMat(:,i)*(xMat(:,i)'*sigma_stimate*xMat(:,i)+tau^2)^-1*(r(i)-xMat(:,i)'*w_stimate);
    if i==n/2+1
       sigmaMat{1,i}(2,2)=sigma0;
    end
    sigma1=[sigma1 sigmaMat{1,i}(1,1)];
    sigma2=[sigma2 sigmaMat{1,i}(2,2)];
end
sigma1=[sigma0 sigma1];
sigma2=sigma2(n/2:end);
end
function [w,sigmaMat,beta]=KalmanFilter2(xMat,n,r,stimuliNum,tau,processNoise,sigma0,gamma)
w=zeros(stimuliNum,n);
sigmaMat=zeros(1,n);
sigmaMat(1)=sigma0;
for i=2:n
    sigma_stimate=sigmaMat(i-1)+processNoise;
    w_stimate=w(i-1);
    sigmaMat(i)=sigma_stimate-sigma_stimate*xMat(i)*(xMat(i)'*sigma_stimate*xMat(i)+tau^2)^-1*xMat(i)'*sigma_stimate;
    w(i)=w_stimate+sigma_stimate*xMat(i)*(xMat(i)'*sigma_stimate*xMat(i)+tau^2)^-1*(r(i)-xMat(i)'*w_stimate);
    beta(i) = (r(i) - xMat(i)*w(i))^2 / (xMat(i)*sigmaMat(i) + tau^2);
    if beta(i) > gamma
            sigmaMat(i) = 40;
    end
end  
end


