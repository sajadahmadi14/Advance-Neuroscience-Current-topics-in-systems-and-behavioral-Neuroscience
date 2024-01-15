%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_4


%% Part 1_a , Weltch
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs=200;

fLog=zeros(48,128);
powerLog=zeros(48,128);
for i=1:length(chan)
    [pxx(i,:),f(i,:)] = pwelch(mean(chan(i).lfp(:,Intersect_Clean_Trials),2),[],[0],[],fs)
    fLog(i,:) = log10(f(i,2:end));
    powerLog(i,:)=log10(pxx(i,2:end));
    coeff(i,:) = polyfit(fLog(i,:),powerLog(i,:),1);
    backGpower(i,:) = coeff(i,1)*fLog(i,:) + coeff(i,2);
    backGpower2(i,:) = coeff(i,1)*log10(f(i,:)) + coeff(i,2);
    NormalizedPower(i,:)=powerLog(i,:)-backGpower(i,:);
end
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=1:16
    subplot(4,4,i)
    plot(fLog(i,:),powerLog(i,:),'k')
    hold on
    plot(fLog(i,:),backGpower(i,:),'r')
    xlim([0,2])
    ylim([-5.1 5.1])
    xlabel('Frequency (Log Scale)','interpreter','Latex')
    ylabel('Power (Log Scale)','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    end
print -depsc fig1.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=17:32
    subplot(4,4,i-16)
    plot(fLog(i,:),powerLog(i,:),'k')
    hold on
    plot(fLog(i,:),backGpower(i,:),'r')
    xlim([0,2])
    ylim([-5.1 5.1])
    xlabel('Frequency (Log Scale)','interpreter','Latex')
    ylabel('Power (Log Scale)','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    end
print -depsc fig2.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=33:48
    subplot(4,4,i-32)
    plot(fLog(i,:),powerLog(i,:),'k')
    hold on
    plot(fLog(i,:),backGpower(i,:),'r')
    xlim([0,2])
    ylim([-5.1 5.1])
    xlabel('Frequency (Log Scale)','interpreter','Latex')
    ylabel('Power (Log Scale)','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    end
print -depsc fig3.eps
    
    figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=1:16
    subplot(4,4,i)
    plot(f(i,2:end),NormalizedPower(i,:),'k')
    hold on
    xlabel('Frequency (Original Scale)','interpreter','Latex')
    ylabel('Normalized Power ','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    index(i)=find(NormalizedPower(i,13:59)==max(NormalizedPower(i,13:59)));
    freq=f(i,:);
    dominantfreq(i)=freq(index(i)+13);
    end
print -depsc fig4.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=17:32
    subplot(4,4,i-16)
    plot(f(i,2:end),NormalizedPower(i,:),'k')
    hold on
    xlabel('Frequency (Original Scale)','interpreter','Latex')
    ylabel('Normalized Power ','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    index(i)=find(NormalizedPower(i,13:59)==max(NormalizedPower(i,13:59)));
    freq=f(i,:);
    dominantfreq(i)=freq(index(i)+13);
    end
print -depsc fig5.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=33:48
    subplot(4,4,i-32)
    plot(f(i,2:end),NormalizedPower(i,:),'k')
    hold on
    xlabel('Frequency (Original Scale)','interpreter','Latex')
    ylabel('Normalized Power','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    index(i)=find(NormalizedPower(i,13:59)==max(NormalizedPower(i,13:59)));
    freq=f(i,:);
    dominantfreq(i)=freq(index(i)+13);
    end
print -depsc fig6.eps

dominantMat=ChannelPosition;
for k=1:48
    dominantMat(find(ChannelPosition==k))=dominantfreq(k);
end
figure;
colormap()
imagesc([1:10],[1:5] ,dominantMat)
set(gca,'YDir','reverse')
c=colorbar
c.Label.String = 'Frequency';
title('Dominant frequencies (Welch)','interpreter', 'latex')
print -depsc ColorCode1.eps

save('dominantfreq.mat','dominantfreq')
    
%% Part 1_a , Multitaper
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs=200;

fLog=zeros(48,320);
powerLog=zeros(48,320);
for i=1:length(chan)
    [pxx(i,:),f(i,:)] = pmtm(mean(chan(i).lfp(:,Intersect_Clean_Trials),2),1.25,length(mean(chan(i).lfp(:,Intersect_Clean_Trials),2)),fs);
    fLog(i,:) = log10(f(i,2:end));
    powerLog(i,:)=log10(pxx(i,2:end));
    coeff(i,:) = polyfit(fLog(i,:),powerLog(i,:),1);
    backGpower(i,:) = coeff(i,1)*fLog(i,:) + coeff(i,2);
    backGpower2(i,:) = coeff(i,1)*log10(f(i,:)) + coeff(i,2);
    NormalizedPower(i,:)=powerLog(i,:)-backGpower(i,:);
end

figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=1:16
    subplot(4,4,i)
    plot(fLog(i,:),powerLog(i,:),'k')
    hold on
    plot(fLog(i,:),backGpower(i,:),'r')
    xlim([0,2])
    ylim([-5.1 5.1])
    xlabel('Frequency (Log Scale)','interpreter','Latex')
    ylabel('Power (Log Scale)','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    end
print -depsc fig7.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=17:32
    subplot(4,4,i-16)
    plot(fLog(i,:),powerLog(i,:),'k')
    hold on
    plot(fLog(i,:),backGpower(i,:),'r')
    xlim([0,2])
    ylim([-5.1 5.1])
    xlabel('Frequency (Log Scale)','interpreter','Latex')
    ylabel('Power (Log Scale)','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    end
print -depsc fig8.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=33:48
    subplot(4,4,i-32)
    plot(fLog(i,:),powerLog(i,:),'k')
    hold on
    plot(fLog(i,:),backGpower(i,:),'r')
    xlim([0,2])
    ylim([-5.1 5.1])
    xlabel('Frequency (Log Scale)','interpreter','Latex')
    ylabel('Power (Log Scale)','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    end
print -depsc fig9.eps
    
    figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=1:16
    subplot(4,4,i)
    plot(f(i,2:end),NormalizedPower(i,:),'k')
    hold on
    xlabel('Frequency (Original Scale)','interpreter','Latex')
    ylabel('Normalized Power ','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    index(i)=find(NormalizedPower(i,34:146)==max(NormalizedPower(i,34:146)));
    freq=f(i,:);
    dominantfreq(i)=freq(index(i)+34);
    end
print -depsc fig10.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=17:32
    subplot(4,4,i-16)
    plot(f(i,2:end),NormalizedPower(i,:),'k')
    hold on
    xlabel('Frequency (Original Scale)','interpreter','Latex')
    ylabel('Normalized Power ','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    index(i)=find(NormalizedPower(i,34:146)==max(NormalizedPower(i,34:146)));
    freq=f(i,:);
    dominantfreq(i)=freq(index(i)+34);
    end
print -depsc fig11.eps
figure;
set(gcf,'units','points','position',[0,0,1200,600])
    for i=33:48
    subplot(4,4,i-32)
    plot(f(i,2:end),NormalizedPower(i,:),'k')
    hold on
    xlabel('Frequency (Original Scale)','interpreter','Latex')
    ylabel('Normalized Power','interpreter','Latex')
    title(['Channel ',num2str(i)],'interpreter', 'latex')
    index(i)=find(NormalizedPower(i,34:146)==max(NormalizedPower(i,34:146)));
    freq=f(i,:);
    dominantfreq(i)=freq(index(i)+34);
    end
print -depsc fig12.eps

dominantMat=ChannelPosition;
for k=1:48
    dominantMat(find(ChannelPosition==k))=dominantfreq(k);
end
figure;
colormap()
imagesc([1:10],[1:5] ,dominantMat)
set(gca,'YDir','reverse')
c=colorbar
c.Label.String = 'Frequency';
title('Dominant frequencies (Multitaper)','interpreter', 'latex')

print -depsc ColorCode2.eps


%% part 1_c , Welch
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs=200;

tVec=1:25:641;
for i=1:length(chan)
for j=1:length(tVec)-1
    [pxx(i,j,:),f] = pwelch(mean(chan(i).lfp(tVec(j):tVec(j+1),Intersect_Clean_Trials),2),[],[],[26],fs)
end
end
pxxAvg=mean(pxx,1);
for k=1:length(f)
pxxAvgNew(k,:)=pxxAvg(:,:,k)
end
figure;
colormap()
imagesc(Time,f,pxxAvgNew)
set(gca,'YDir','normal')
c=colorbar
c.Label.String = 'Power';
title('Welch','interpreter', 'latex')
xlabel('Time(s)','interpreter','Latex')
ylabel('Frequency (Hz) ','interpreter','Latex')
print -depsc PowSpec1.eps

%% part 1_c , Multitaper
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs=200;

tVec=1:25:641;
for i=1:length(chan)
for j=1:length(tVec)-1
    [pxx(i,j,:),f] = pmtm(mean(chan(i).lfp(tVec(j):tVec(j+1),Intersect_Clean_Trials),2),[],length(mean(chan(i).lfp(tVec(j):tVec(j+1),Intersect_Clean_Trials),2)),fs);
end
end

pxxAvg=mean(pxx,1);
for k=1:length(f)
pxxAvgNew(k,:)=pxxAvg(:,:,k)
end
figure;
colormap()
imagesc(Time,f,pxxAvgNew)
set(gca,'YDir','normal')
c=colorbar
c.Label.String = 'Power';
title('Multitaper','interpreter', 'latex')
xlabel('Time(s)','interpreter','Latex')
ylabel('Frequency (Hz) ','interpreter','Latex')
print -depsc PowSpec2.eps
  
%% part 2_a
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs=200;

load('dominantfreq.mat')
for i=1:length(chan)
    for nTrials=1:length(Intersect_Clean_Trials)
        [b,a] = butter(2,[dominantfreq(i)-1 dominantfreq(i)+1]/100,'bandpass');
        filteredLFP(i).lfp(:,nTrials)=filtfilt(b,a,chan(i).lfp(:,nTrials));
    end
end

for k=1:length(chan)
for j=1:length(Intersect_Clean_Trials)
        cleanLFP(k).lfp(:,j)=chan(k).lfp(:,Intersect_Clean_Trials(j));
end
end

 % calculate instantaneous phase
 for i=1:length(chan)
    for nTrials=1:length(Intersect_Clean_Trials)
        HilbertSignal=hilbert(filteredLFP(i).lfp(:,nTrials));
        phi(i).phase(:,nTrials)=angle(filteredLFP(i).lfp(:,nTrials)+1j*HilbertSignal);
    end
 end
 save('Phi.mat','phi');
 % design a demo
 t=-1.2:0.005:2;
 NumTrial=randi(490);
 NumTrial=89;
 for q=1:length(t)
     demoMat(q).cosPhi=ChannelPosition;
 for k=1:48
     demoMat(q).cosPhi(find(ChannelPosition==k))=cos(phi(k).phase(q,NumTrial));
 end
 end
 
 slot=[312:319];
 figure;
 set(gcf,'units','points','position',[0,0,1200,300])
 for n=1:length(slot)
     subplot(2,4,n)
     colormap(jet)
     imagesc(demoMat(slot(n)).cosPhi)
     set(gca,'YDir','normal')
     c=colorbar
     c.Label.String = 'cos(\phi)';
     title(['t=',num2str(t(slot(n))),' s'],'interpreter', 'latex')
 end
 %print -depsc trav.eps

 %% part 2_d : PGD, Speed , Direction of Propagation
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
load('Phi.mat')
fs=200;
t=-1.2:0.005:2;
temp1 = [];
temp2 = [];
PGD=[];
Speed=[];
pgd=[];
speed=[];
DirectionPropag=[];
 for NumTrial=1:490
 for q=1:length(t)
     demoMat(q).cosPhi=ChannelPosition;
 for k=1:48
     demoMat(q).cosPhi(find(ChannelPosition==k))=cos(phi(k).phase(q,NumTrial));
 end
 end
 
 for a=1:641
     [Fx(a).grad,Fy(a).grad]=gradient(demoMat(a).cosPhi(:,2:end));
 end 
 for a=1:640
     speedNumerator(a)=sqrt((mean(mean(Fx(a+1).grad))-mean(mean(Fx(a).grad)))^2+(mean(mean(Fy(a+1).grad))-mean(mean(Fy(a).grad)))^2)*fs;
     for i=1:size(Fx(a).grad,1)
         for j=1:size(Fx(a).grad,2)
             denominator(i,j)=norm([Fx(a).grad(i,j),Fy(a).grad(i,j)]);
         end
     end
     pgd(a)=norm([mean(mean(Fx(a).grad)),mean(mean(Fy(a).grad))])/mean(mean(denominator));
     speed(a)=speedNumerator(a)/mean(mean(denominator));
 end
 PGD=[PGD,pgd];
 Speed=[Speed,speed];
      for m=1:640
          temp1 = [temp1, reshape(Fx(m).grad,1,size(Fx(a).grad,1)*size(Fx(a).grad,2))];
          temp2 = [temp2, reshape(Fy(m).grad,1,size(Fx(a).grad,1)*size(Fx(a).grad,2))];
      end
 end
 
PGDMat=reshape(PGD,[640,313]);
SpeedMat=reshape(Speed,[640,313]);

PGDavg=mean(PGDMat,2);
Speedavg=mean(SpeedMat,2);

for i = 1:length(temp1)
     x=temp1(1,i); y=temp2(1,i);
     if (x>=0) 
         DirectionPropag(i) = atand(y/x);
     elseif (x<0 && y>=0)
         DirectionPropag(i) = atand(y/x) + 180;
     elseif (x<0 && y<0)
         DirectionPropag(i) = atand(y/x) - 180;
     end
     if (DirectionPropag(i) <= -180)
         DirectionPropag(i) = -DirectionPropag(i);
     end
end

xPGD=min(PGD):0.003:max(PGD);
hPGD = hist(PGD , length(xPGD));
figure;
bar(x,hPGD,'b','linewidth',8);
xlabel('PGD','interpreter', 'latex');
ylabel('Count','interpreter', 'latex');
title('Histogram of PGD during the instruction epoch of all trials for all times','interpreter', 'latex')
print -depsc PGDhis.eps
    
xSpeed=min(Speed):0.8:max(Speed);
hSpeed = hist(Speed , length(xSpeed));
figure;
bar(xSpeed,hSpeed,'b','linewidth',8);
xlim([0,150])
xlabel('Speed (cm/s)','interpreter', 'latex');
ylabel('Count','interpreter', 'latex');
title('Histogram of Speed during the instruction epoch of all trials for all times','interpreter', 'latex')
print -depsc Speedhis.eps

figure;
set(gcf,'units','points','position',[0,0,1200,300])
plot([-1.2:0.005:2-0.005],PGDavg,'b');
xlabel('t(s)','interpreter', 'latex');
ylabel('PGD','interpreter', 'latex');
title('PGD averaged across all trials as a function of time','interpreter', 'latex')
xlim([-1.2,2])
print -depsc pgdvstime.eps

figure;
set(gcf,'units','points','position',[0,0,1200,300])
plot([-1.2:0.005:2-0.005],Speedavg,'b');
xlabel('t(s)','interpreter', 'latex');
ylabel('Speed (cm/s)','interpreter', 'latex');
title('Speed averaged across all trials as a function of time','interpreter', 'latex')
xlim([-1.2,2])
print -depsc speedvstime.eps

figure;
set(gcf,'units','points','position',[0,0,800,400])
direction = -180:1:180; 
hdirec = hist(DirectionPropag,length(direction));
prefDirection=direction(hdirec==max(hdirec));
bar(direction,hdirec,'b')
title('Histogram of Direction of Propagation during the instruction epoch of all trials for all times','interpreter','latex')
xlabel('Propagation direction(deg)','interpreter','latex')
ylabel('Direction count','interpreter','latex')
print -depsc progDirec.eps



