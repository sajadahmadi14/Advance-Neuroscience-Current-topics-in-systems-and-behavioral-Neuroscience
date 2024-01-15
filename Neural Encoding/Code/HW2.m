

%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_2

%%
clear;
clc;
close all;

load ('UnitsData.mat')
%% Step 1
% Psth plot for 16 randomly selected unit (Figure1)
Trials=1:192;
window_width=0.1; 
%unit=randi([1,481],1,20);
unit=[10 17 245 47 400 178 94 347 4 80 214 420 64 472 32 210 380 297 100 409]
%unit=1:481;
figure;
set(gcf,'units','points','position',[0,0,1200,600])
for j=1:length(Trials)
[Psth,tvec]=PsthFunc(Unit,window_width,unit(i),Trials);
subplot(5,4,i)
plot(tvec,Psth,'color','b')
xlim([min(tvec) max(tvec)])
ylim([min(Psth)-2 max(Psth)+2])
xlabel('Time (s)','interpreter', 'latex')
ylabel('FR(Hz)','interpreter', 'latex');
title(['PSTH for Unit ',num2str(unit(i))],'fontsize',8,'interpreter', 'latex')
end
%print -depsc fig1.eps

% Calculate the average psth for different conditions(1 to 6) (Figure2)
figure;
set(gcf,'units','points','position',[0,0,1200,600])
color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4940 0.1840 0.5560],[0 0 1],[0.75 0 0],[0.25 0 0]};
for k=1:length(unit)
window_width=0.1;
for i=1:6
        Trials=Unit(unit(k)).Cnd(i).TrialIdx;
        [psth(i,:),tvec(i,:)]=PsthFunc(Unit,window_width,unit(k),Trials);
        subplot(5,4,k)
        plot(tvec(i,:),psth(i,:),'color',color{1,i},'LineWidth',1.4);
        xlim([min(tvec(i,:)) max(tvec(i,:))])
        xlabel('Time (s)','interpreter', 'latex')
        ylabel('FR(Hz)','interpreter', 'latex');
        title(['Average PSTH for Unit ',num2str(unit(k))],'fontsize',8,'interpreter', 'latex')
        hold on;
end
end
h=get(gca,'Children');
fig = gcf;
fig.Position(3) = fig.Position(3) + 300;
Lgnd=legend("Condition1","Condition2","Condition3","Condition4","Condition5","Condition6",'interpreter', 'latex')
Lgnd.Position(1) = 0.02;
Lgnd.Position(2) = 0.6;
%print -depsc fig2.eps

% Calculate the average psth for reward conditions(3,6,9) (Figure3)
unit=[10 17 245 47 400 178 94 347 4 80 214 420 64 472 32 210 380 297 100 409]
figure;
set(gcf,'units','points','position',[0,0,1200,600])
for k=1:length(unit)
window_width=0.1;
for i=1:6
          Trials=Unit(unit(k)).Cnd(i).TrialIdx;
         [psth(i,:),tvec(i,:)]=PsthFunc(Unit,window_width,unit(k),Trials);
end
psth3=(psth(1,:)+psth(2,:))/2;
psth6=(psth(3,:)+psth(4,:))/2;
psth9=(psth(5,:)+psth(6,:))/2;

subplot(5,4,k)
plot(tvec,psth3,'color',[0.6350 0.0780 0.1840],'LineWidth',1.4);
hold on
plot(tvec,psth6,'color',[0 0.4470 0.7410],'LineWidth',1.5);
plot(tvec,psth9,'color',[0.4940 0.1840 0.5560],'LineWidth',1.5);
hold off
xlim([min(tvec(i,:)) max(tvec(i,:))]);
xlabel('Time (s)','interpreter', 'latex')
ylabel('FR(Hz)','interpreter', 'latex');
title(['Average PSTH for Unit ',num2str(unit(k))],'fontsize',8,'interpreter', 'latex')
end

h=get(gca,'Children');
fig = gcf;
fig.Position(3) = fig.Position(3) + 300;
Lgnd=legend(h([1 7 13]),"reward=3","reward=6","reward=9",'interpreter', 'latex');
Lgnd.Position(1) = 0.02;
Lgnd.Position(2) = 0.6;
%print -depsc fig3.eps

% Calculate the average psth for cue conditions(-1,1) (Figure4)
unit=[10 17 245 47 400 178 94 347 4 80 214 420 64 472 32 210 380 297 100 409];
figure;
set(gcf,'units','points','position',[0,0,1200,600])
for k=1:length(unit)
window_width=0.1;
for i=1:6
          Trials=Unit(unit(k)).Cnd(i).TrialIdx;
         [psth(i,:),tvec(i,:)]=PsthFunc(Unit,window_width,unit(k),Trials);
end
psth1=(psth(2,:)+psth(4,:)+psth(6,:))/3;
psth_1=(psth(1,:)+psth(3,:)+psth(5,:))/3;

subplot(5,4,k)
plot(tvec,psth1,'color',[0.6350 0.0780 0.1840],'LineWidth',1.4);
hold on
plot(tvec,psth_1,'color',[0 0.4470 0.7410],'LineWidth',1.5);
hold off
xlim([min(tvec(i,:)) max(tvec(i,:))]);
xlabel('Time (s)','interpreter', 'latex')
ylabel('FR(Hz)','interpreter', 'latex');
title(['Average PSTH for Unit ',num2str(unit(k))],'fontsize',8,'interpreter', 'latex')
end

h=get(gca,'Children');
fig = gcf;
fig.Position(3) = fig.Position(3) + 300;
Lgnd=legend(h([1 7 ]),"cue=1","cue=-1",'interpreter', 'latex');
Lgnd.Position(1) = 0.02;
Lgnd.Position(2) = 0.6;
print -depsc fig4.eps
%% Part 2
%% P-value , Reward=3,6,9 , WindowSize=50ms
for i=1:length(unit)    
[Psth(i).trials,tvec]=PsthFix(Unit,window_width,unit(i),Trials);
end
for i=1:481
    cnd3(i).trials=[Unit(i).Cnd(1).TrialIdx;Unit(i).Cnd(2).TrialIdx]
    cnd6(i).trials=[Unit(i).Cnd(3).TrialIdx;Unit(i).Cnd(4).TrialIdx]
    cnd9(i).trials=[Unit(i).Cnd(5).TrialIdx;Unit(i).Cnd(6).TrialIdx]
end

CndTrials=zeros(481,192);
for i=1:481 
for j=1:192
    CndTrials(i,cnd3(i).trials)=3
    CndTrials(i,cnd6(i).trials)=6
    CndTrials(i,cnd9(i).trials)=9
end
end

for i=1:481
mdl = fitglm(Psth(i).trials,CndTrials(i,:))
pValueUnit(i,1)=coefTest(mdl);
pValueAll(i,:)=mdl.Coefficients.pValue(2:65);
end

unit=[5,12,30,48]
figure('NumberTitle', 'off', 'Name', 'P-Value for units (Reward=3 , Reward=6 , Reward=9)')
set(gcf,'units','points','position',[0,0,1200,600])
imagesc(pValueAll(unit,:))
c=colorbar
c.Label.String = 'P-Value','interpreter', 'latex';
xlabel('Window (each window is 50 ms)','interpreter', 'latex')
ylabel('Unit','interpreter', 'latex');
title('P-Value for units (5,12,30,40) ','fontsize',8,'interpreter', 'latex')
%print -depsc fig5_b.eps

figure('NumberTitle', 'off', 'Name', 'P-Value for units (Reward=3 , Reward=6 , Reward=9)')
set(gcf,'units','points','position',[0,0,1200,600])
m=[1 61 121 181 241 301 361 421 481]
n=[60 120 180 240 300 360 420 480]
for i=1:length(n)
    subplot(4,2,i)
imagesc(pValueAll(m(i):n(i),:))
c=colorbar
c.Label.String = 'P-Value','interpreter', 'latex';
xlabel('Window (each window is 50 ms)','interpreter', 'latex')
ylabel('Unit','interpreter', 'latex');
title(['P-Value for units ',num2str(m(i)),'to',num2str(n(i))],'fontsize',8,'interpreter', 'latex')
end
%print -depsc fig5.eps

unit=[5 12 30 48]
figure;
set(gcf,'units','points','position',[0,0,600,350])
for k=1:length(unit)
window_width=0.03;
for i=1:6
          Trials=Unit(unit(k)).Cnd(i).TrialIdx;
         [psth(i,:),tvec(i,:)]=PsthFunc(Unit,window_width,unit(k),Trials);
end
psth3=(psth(1,:)+psth(2,:))/2;
psth6=(psth(3,:)+psth(4,:))/2;
psth9=(psth(5,:)+psth(6,:))/2;
subplot(2,2,k)
plot(tvec,psth3,'color',[0.6350 0.0780 0.1840],'LineWidth',1.4);
hold on
plot(tvec,psth6,'color',[0 0.4470 0.7410],'LineWidth',1.5);
plot(tvec,psth9,'color',[0.4940 0.1840 0.5560],'LineWidth',1.5);
hold off
xlim([min(tvec(i,:)) max(tvec(i,:))]);
xlabel('Time (s)','interpreter', 'latex')
ylabel('FR(Hz)','interpreter', 'latex');
title(['Average PSTH for Unit ',num2str(unit(k))],'fontsize',8,'interpreter', 'latex')
end

h=get(gca,'Children');
fig = gcf;
fig.Position(3) = fig.Position(3) + 300;
Lgnd=legend(h([1 7 13]),"reward=3","reward=6","reward=9",'interpreter', 'latex');
Lgnd.Position(1) = 0.0;
Lgnd.Position(2) = 0.6;
%print -depsc fig5_c.eps

%% P-value , Direction=1,-1 , WindowSize=50ms

for i=1:481
    cnd1(i).trials=[Unit(i).Cnd(2).TrialIdx;Unit(i).Cnd(4).TrialIdx;Unit(i).Cnd(6).TrialIdx]
    cnd_1(i).trials=[Unit(i).Cnd(1).TrialIdx;Unit(i).Cnd(3).TrialIdx;Unit(i).Cnd(5).TrialIdx]
end

CndTrials=zeros(481,192);
for i=1:481 
    CndTrials(i,cnd1(i).trials)=1
    CndTrials(i,cnd_1(i).trials)=-1
end
for i=1:481
mdl = fitglm(Psth(i).trials,CndTrials(i,:))
pValueUnit(i,1)=coefTest(mdl);
pValueAll(i,:)=mdl.Coefficients.pValue(2:65);
end

figure('NumberTitle', 'off', 'Name', 'P-Value for units (Condition=1,-1)')
set(gcf,'units','points','position',[0,0,1200,600])
m=[1 61 121 181 241 301 361 421 481]
n=[60 120 180 240 300 360 420 480]
for i=1:length(n)
    subplot(4,2,i)
imagesc(pValueAll(m(i):n(i),:))
c=colorbar
c.Label.String = 'P-Value','interpreter', 'latex';
xlabel('Window (each window is 50 ms)','interpreter', 'latex')
ylabel('Unit','interpreter', 'latex');
title(['P-Value for units ',num2str(m(i)),'to',num2str(n(i))],'fontsize',8,'interpreter', 'latex')
end
%print -depsc fig7.eps

figure('NumberTitle', 'off', 'Name', 'P-Value for units (Condition=1,-1)')
unit=[13 16 27 40]
set(gcf,'units','points','position',[0,0,1200,600])
imagesc(pValueAll(unit,:))
c=colorbar
c.Label.String = 'P-Value','interpreter', 'latex';
xlabel('Window (each window is 50 ms)','interpreter', 'latex')
ylabel('Unit','interpreter', 'latex');
title('P-Value for units(13,16,27,40) ','fontsize',8,'interpreter', 'latex')
%print -depsc fig7_b.eps

unit=[13 16 27 40];
figure;
set(gcf,'units','points','position',[0,0,500,350])
for k=1:length(unit)
window_width=0.03;
for i=1:6
          Trials=Unit(unit(k)).Cnd(i).TrialIdx;
         [psth(i,:),tvec(i,:)]=PsthFunc(Unit,window_width,unit(k),Trials);
end
psth1=(psth(2,:)+psth(4,:)+psth(6,:))/3;
psth_1=(psth(1,:)+psth(3,:)+psth(5,:))/3;
subplot(2,2,k)
plot(tvec,psth1,'color',[0.6350 0.0780 0.1840],'LineWidth',1.4);
hold on
plot(tvec,psth_1,'color',[0 0.4470 0.7410],'LineWidth',1.5);
hold off
xlim([min(tvec(i,:)) max(tvec(i,:))]);
xlabel('Time (s)','interpreter', 'latex')
ylabel('FR(Hz)','interpreter', 'latex');
title(['Average PSTH for Unit ',num2str(unit(k))],'fontsize',8,'interpreter', 'latex')
end

h=get(gca,'Children');
fig = gcf;
fig.Position(3) = fig.Position(3) + 300;
Lgnd=legend(h([1 7 ]),"cue=1","cue=-1",'interpreter', 'latex');
Lgnd.Position(1) = 0.0;
Lgnd.Position(2) = 0.6;
%print -depsc fig8.eps

%% Part 3
window_width=0.03
       for cnd=1:6
            for nUnit=1:481
            ConditionTrials=Unit(nUnit).Cnd(cnd).TrialIdx  
            [Cnd(cnd).Psth(nUnit,:)]=PsthFunc(Unit,window_width,nUnit,ConditionTrials);
            end
       end
  % PCA (Condition 1 to 6)
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(Cnd(1).Psth(:,1200:3199)')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(Cnd(2).Psth(:,1200:3199)')
  [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(Cnd(3).Psth(:,1200:3199)')
  [coeff4,score4,latent4,tsquared4,explained4,mu4] = pca(Cnd(4).Psth(:,1200:3199)')
  [coeff5,score5,latent5,tsquared5,explained5,mu5] = pca(Cnd(5).Psth(:,1200:3199)')
  [coeff6,score6,latent6,tsquared6,explained6,mu6] = pca(Cnd(6).Psth(:,1200:3199)')

 figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],2,1,3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],2,1,3) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],2,1,3) 
 hold on 
 PCAplot(score4,[0 0 1],2,1,3) 
 hold on 
 PCAplot(score5,[1 0 0],2,1,3) 
 hold on 
 PCAplot(score6,[0 0 0],2,1,3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([18 15 12 9 6 3 ]),"Condition1","Condition2","Condition3","Condition4","Condition5","Condition6",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig18.eps

 figure;
  PCAplot(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],3) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],3) 
 hold on 
 PCAplot(score4,[0 0 1],3) 
 hold on 
 PCAplot(score5,[1 0 0],3) 
 hold on 
 PCAplot(score6,[0 0 0],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([18 15 12 9 6 3 ]),"Condition1","Condition2","Condition3","Condition4","Condition5","Condition6",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig10.eps

% PCA (Reward 3,6,9)
CND3=(Cnd(1).Psth(:,1200:3199)+Cnd(2).Psth(:,1200:3199))/2;
CND6=(Cnd(3).Psth(:,1200:3199)+Cnd(4).Psth(:,1200:3199))/2;
CND9=(Cnd(5).Psth(:,1200:3199)+Cnd(6).Psth(:,1200:3199))/2;
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(CND3')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(CND6')
  [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(CND9')
   figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],2,1,3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],2,1,3) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],2,1,3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([9 6 3 ]),"Reward=3","Reward=6","Reward=9",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%  print -depsc fig19.eps
 
 figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],3) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([9 6 3 ]),"Reward=3","Reward=6","Reward=9",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig14.eps
  
% CUE 1,-1
CND1=(Cnd(2).Psth(:,1200:3199)+Cnd(4).Psth(:,1200:3199)+Cnd(6).Psth(:,1200:3199))/3;
CND_1=(Cnd(1).Psth(:,1200:3199)+Cnd(3).Psth(:,1200:3199)+Cnd(5).Psth(:,1200:3199))/3;
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(CND1')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(CND_1')
 figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],2,1,3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],2,1,3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([ 6 3 ]),"Cue=1","Cue=-1",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig20.eps

 figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([ 6 3 ]),"Cue=1","Cue=-1",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig17.eps

%% Part 4
%% CFR
startup
load ('UnitsData.mat')
t=-1.2:0.05:2-0.05
window_width=0.05;
unit=[1:481]
for k=1:length(unit)
for i=1:6
          Trials=Unit(unit(k)).Cnd(i).TrialIdx;
         [PSTH(i,k,:)]=PsthFix2(Unit,window_width,unit(k),Trials);
end
end
dataTensor=permute(PSTH,[3 2 1]);
rng('shuffle', 'twister') 
surrogate_type = 'surrogate-TNC';
model_dim = 6;
times_msk = t>0 & t<2; 
[targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures(dataTensor);
numSurrogates = 100;
params = [];
params.readout_mode = 2;         
params.shfl_mode = 3;        
params.fix_mode = 2;        

if strcmp(surrogate_type, 'surrogate-T')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = [];
    params.margCov{3} = [];
    params.meanTensor = M.T;
elseif strcmp(surrogate_type, 'surrogate-TN')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = [];
    params.meanTensor = M.TN;
elseif strcmp(surrogate_type, 'surrogate-TNC')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.TNC; 
else
    error('please specify a correct surrogate type') 
end
for i = 1:numSurrogates
    [surrTensor] = sampleCFR(dataTensor, params)     
end
%save('suurTensorCFR.mat','surrTensor')

load('suurTensorCFR.mat')
load('Cnd.mat')
  psth=permute(surrTensor,[3 2 1])
  for i=1:6
     psthTensor(i).cnd=reshape(psth(i,:,:),481,64) 
  end
  % PCA (Condition 1 to 6)
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(psthTensor(1).cnd(:,24:63)')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(psthTensor(2).cnd(:,24:63)')
  [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(psthTensor(3).cnd(:,24:63)')
  [coeff4,score4,latent4,tsquared4,explained4,mu4] = pca(psthTensor(4).cnd(:,24:63)')
  [coeff5,score5,latent5,tsquared5,explained5,mu5] = pca(psthTensor(5).cnd(:,24:63)')
  [coeff6,score6,latent6,tsquared6,explained6,mu6] = pca(psthTensor(6).cnd(:,24:63)')

 figure;
 PCAplotshuffle(score1,[0.4940 0.1840 0.5560],2,3,2) 
 hold on 
 PCAplotshuffle(score2,[0.4660 0.6740 0.1880],2,3,2) 
 hold on 
 PCAplotshuffle(score3,[0.6350 0.0780 0.1840],2,3,2) 
 hold on 
 PCAplotshuffle(score4,[0 0 1],2,3,2) 
 hold on 
 PCAplotshuffle(score5,[1 0 0],2,3,2) 
 hold on 
 PCAplotshuffle(score6,[0 0 0],2,3,2) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([18 15 12 9 6 3 ]),"Condition1","Condition2","Condition3","Condition4","Condition5","Condition6",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig23.eps

 figure;
  PCAplotshuffle(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplotshuffle(score2,[0.4660 0.6740 0.1880],3) 
 hold on 
 PCAplotshuffle(score3,[0.6350 0.0780 0.1840],3) 
 hold on 
 PCAplotshuffle(score4,[0 0 1],3) 
 hold on 
 PCAplotshuffle(score5,[1 0 0],3) 
 hold on 
 PCAplotshuffle(score6,[0 0 0],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([18 15 12 9 6 3 ]),"Condition1","Condition2","Condition3","Condition4","Condition5","Condition6",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig24.eps

% PCA (Reward 3,6,9)
CND3=(psthTensor(1).cnd(:,24:63)+psthTensor(2).cnd(:,24:63))/2;
CND6=(psthTensor(3).cnd(:,24:63)+psthTensor(4).cnd(:,24:63))/2;
CND9=(psthTensor(5).cnd(:,24:63)+psthTensor(6).cnd(:,24:63))/2;
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(CND3')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(CND6')
  [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(CND9')
   figure;
 PCAplotshuffle(score1,[0.4940 0.1840 0.5560],2,3,2) 
 hold on 
 PCAplotshuffle(score2,[0.4660 0.6740 0.1880],2,3,2) 
 hold on
 PCAplotshuffle(score3,[0.6350 0.0780 0.1840],2,3,2) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([9 6 3 ]),"Reward=3","Reward=6","Reward=9",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
 % print -depsc fig27.eps
 
   figure;
 PCAplotshuffle(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplotshuffle(score2,[0.4660 0.6740 0.1880],3) 
 hold on 
 PCAplotshuffle(score3,[0.6350 0.0780 0.1840],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([9 6 3 ]),"Reward=3","Reward=6","Reward=9",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig28.eps
  
% CUE 1,-1
CND1=(psthTensor(2).cnd(:,24:63)+psthTensor(4).cnd(:,24:63)+psthTensor(6).cnd(:,24:63))/3;
CND_1=(psthTensor(1).cnd(:,24:63)+psthTensor(3).cnd(:,24:63)+psthTensor(5).cnd(:,24:63))/3;
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(CND1')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(CND_1')
   figure;
 PCAplotshuffle(score1,[0.4940 0.1840 0.5560],2,3,2) 
 hold on 
 PCAplotshuffle(score2,[0.4660 0.6740 0.1880],2,3,2) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([ 6 3 ]),"Cue=1","Cue=-1",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig31.eps
 figure;
 PCAplotshuffle(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplotshuffle(score2,[0.4660 0.6740 0.1880],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([ 6 3 ]),"Cue=1","Cue=-1",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig32.eps


%% TME

startup
load ('UnitsData.mat')
t=-1.2:0.001:2+0.001
window_width=0.05;
unit=[1:481]
for k=1:length(unit)
for i=1:6
          Trials=Unit(unit(k)).Cnd(i).TrialIdx;
         [PSTH(i,k,:)]=PsthFix2(Unit,window_width,unit(k),Trials);
end
end
dataTensor=permute(PSTH,[3 2 1]);
rng('shuffle', 'twister') % randomize the seed
surrogate_type = 'surrogate-TNC';
model_dim = 6;
times_msk = t>0 & t<2; 
[targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures(dataTensor);
numSurrogates = 100;
params = [];   
if strcmp(surrogate_type, 'surrogate-T')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = [];
    params.margCov{3} = [];
    params.meanTensor = M.T;
elseif strcmp(surrogate_type, 'surrogate-TN')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = [];
    params.meanTensor = M.TN;
elseif strcmp(surrogate_type, 'surrogate-TNC')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.TNC; 
else
    error('please specify a correct surrogate type') 
end
maxEntropy=fitMaxEntropy(params);
for i = 1:numSurrogates
    [surrTensor] = sampleTME(maxEntropy)      
end
save('surrTensorTME.mat','surrTensor')

load('surrTensorTME.mat')
load('Cnd.mat')
  psth=permute(surrTensor,[3 2 1])
  for i=1:6
     psthTensor(i).cnd=reshape(psth(i,:,:),481,3202) 
  end
   PCA (Condition 1 to 6)
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(psthTensor(1).cnd(:,1200:3199)')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(psthTensor(2).cnd(:,1200:3199)')
  [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(psthTensor(3).cnd(:,1200:3199)')
  [coeff4,score4,latent4,tsquared4,explained4,mu4] = pca(psthTensor(4).cnd(:,1200:3199)')
  [coeff5,score5,latent5,tsquared5,explained5,mu5] = pca(psthTensor(5).cnd(:,1200:3199)')
  [coeff6,score6,latent6,tsquared6,explained6,mu6] = pca(psthTensor(6).cnd(:,1200:3199)')

 figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],2,3,2) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],2,3,2) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],2,3,2) 
 hold on 
 PCAplot(score4,[0 0 1],2,3,2) 
 hold on 
 PCAplot(score5,[1 0 0],2,3,2) 
 hold on 
 PCAplot(score6,[0 0 0],2,3,2) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([18 15 12 9 6 3 ]),"Condition1","Condition2","Condition3","Condition4","Condition5","Condition6",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig35.eps

 figure;
  PCAplot(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],3) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],3) 
 hold on 
 PCAplot(score4,[0 0 1],3) 
 hold on 
 PCAplot(score5,[1 0 0],3) 
 hold on 
 PCAplot(score6,[0 0 0],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([18 15 12 9 6 3 ]),"Condition1","Condition2","Condition3","Condition4","Condition5","Condition6",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig36.eps

 
% PCA (Reward 3,6,9)
CND3=(psthTensor(1).cnd(:,1200:3199)+psthTensor(2).cnd(:,1200:3199))/2;
CND6=(psthTensor(3).cnd(:,1200:3199)+psthTensor(4).cnd(:,1200:3199))/2;
CND9=(psthTensor(5).cnd(:,1200:3199)+psthTensor(6).cnd(:,1200:3199))/2;
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(CND3')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(CND6')
  [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(CND9')
  figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],2,3,2) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],2,3,2) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],2,3,2) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([9 6 3 ]),"Reward=3","Reward=6","Reward=9",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
 %print -depsc fig39.eps
 
   figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],3) 
 hold on 
 PCAplot(score3,[0.6350 0.0780 0.1840],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([9 6 3 ]),"Reward=3","Reward=6","Reward=9",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig40.eps
  
% CUE 1,-1
CND1=(psthTensor(2).cnd(:,1200:3199)+psthTensor(4).cnd(:,1200:3199)+psthTensor(6).cnd(:,1200:3199))/3;
CND_1=(psthTensor(1).cnd(:,1200:3199)+psthTensor(3).cnd(:,1200:3199)+psthTensor(5).cnd(:,1200:3199))/3;
  [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(CND1')
  [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(CND_1')
   figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],2,3,2) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],2,3,2) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([ 6 3 ]),"Cue=1","Cue=-1",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig43.eps
  figure;
 PCAplot(score1,[0.4940 0.1840 0.5560],3) 
 hold on 
 PCAplot(score2,[0.4660 0.6740 0.1880],3) 
 legend
 h=get(gca,'Children');
 fig = gcf;
 fig.Position(3) = fig.Position(3) + 300;
 Lgnd=legend(h([ 6 3 ]),"Cue=1","Cue=-1",'interpreter', 'latex');
 Lgnd.Position(1) = 0.75;
 Lgnd.Position(2) = 0.7;
%print -depsc fig44.eps

function [psth,tvec]=PsthFunc(Unit,window_width,unit,Trials)
dt=0.001;
tmax=3.2;
t=0:dt:tmax;
Nt=length(t);
Nwidth = round(window_width/dt);   
two_bw_sq = 2*(window_width)*(window_width); 
tbins = dt*[0:Nt];                    
Nbins = length(tbins);               
psth = zeros(1,Nbins);  
for trial = 1:length(Trials);            
    for i = (floor((Unit(unit).Trls{Trials(trial),1}+1.2)*1000))'     
        prefactor = true;
        for bin = 1:Nbins             
            tdiff = dt*(i-bin+0.5);   
            psth(bin) = psth(bin) + prefactor*exp(-tdiff*tdiff/two_bw_sq);
        end
    end
end
for i = 1:Nbins
    psth(i) = psth(i)/(length(Trials)*sqrt(pi*two_bw_sq) ...
        *(0.5*(erf(tbins(i)/window_width)+erf((tmax-tbins(i))/window_width))));
end
tvec = -1.2:dt:tmax+dt-1.2;                
end

function [psth,tvec]=PsthFix(Unit,window_width,unit,Trials)
dt=0.001;
tmax=3.2;
t=0:dt:tmax;
Nt=length(t);
Nwidth = round(window_width/dt);      
Nwindows1b = round(Nt/Nwidth);         
tvec = dt*[Nwidth/2:Nwidth:Nt];    
psth = zeros(1,Nwindows1b);         
for j=1:length(Trials)
    SpikeTimes=Unit(unit).Trls{Trials(j),1}
for i = 1:Nwindows1b
    psth(j,i) = length(find(-1.2+(i-1)*(Nwidth*dt)<=SpikeTimes & SpikeTimes <= -1.2+(i*(Nwidth*dt))));
end
end
psth = psth/(window_width);  
end

function [psth,tvec]=PsthFix2(Unit,window_width,unit,Trials)
dt=0.001;
tmax=3.2;
t=0:dt:tmax;
Nt=length(t);
Nwidth = round(window_width/dt);       
Nwindows1b = round(Nt/Nwidth);        
tvec = dt*[Nwidth/2:Nwidth:Nt];    
psth = zeros(1,Nwindows1b);       
for j=1:length(Trials)
    SpikeTimes=Unit(unit).Trls{Trials(j),1}
for i = 1:Nwindows1b
    psth(j,i) = length(find(-1.2+(i-1)*(Nwidth*dt)<=SpikeTimes & SpikeTimes <= -1.2+(i*(Nwidth*dt))));
end
end
psth=sum(psth,1);
psth = psth/(length(Trials)*window_width);  
end




function [p,c,r]=PCAplot(score,color,dim,PC1,PC2)
if dim==2
p=plot(score(:,PC1),score(:,PC2),'color',color)
hold on
c=plot(score(1,PC1),score(1,PC2),'o','MarkerFaceColor',color);
hold on
r=plot(score(2000,PC1),score(2000,PC2),'s','MarkerFaceColor',color);
axis square;
xlabel(['Principle axis ',num2str(PC1)],'interpreter', 'latex');
ylabel(['Principle axis ',num2str(PC2)],'interpreter', 'latex');
end
if dim==3
 [x,y,z] = sphere
  p=plot3(score(:,1),score(:,2),score(:,3),'color',color)
  hold on
  c=plot3(score(1,1),score(1,2),score(1,3),'o','MarkerFaceColor',color)
  hold on
  r=plot3(score(2000,1),score(2000,2),score(2000,3),'s','MarkerFaceColor',color) 
  axis square;
  xlabel(['Principle axis ',num2str(1)],'interpreter', 'latex');
  ylabel(['Principle axis ',num2str(2)],'interpreter', 'latex');
  zlabel(['Principle axis ',num2str(3)],'interpreter', 'latex');
end
end

function [p,c,r]=PCAplotshuffle(score,color,dim,PC1,PC2)
if dim==2
p=plot(score(:,PC1),score(:,PC2),'color',color)
hold on
c=plot(score(1,PC1),score(1,PC2),'o','MarkerFaceColor',color);
hold on
r=plot(score(40,PC1),score(40,PC2),'s','MarkerFaceColor',color);
axis square;
xlabel(['Principle axis ',num2str(PC1)],'interpreter', 'latex');
ylabel(['Principle axis ',num2str(PC2)],'interpreter', 'latex');
end
if dim==3
 [x,y,z] = sphere
  p=plot3(score(:,1),score(:,2),score(:,3),'color',color)
  hold on
  c=plot3(score(1,1),score(1,2),score(1,3),'o','MarkerFaceColor',color)
  hold on
  r=plot3(score(40,1),score(40,2),score(40,3),'s','MarkerFaceColor',color) 
  axis square;
  xlabel(['Principle axis ',num2str(1)],'interpreter', 'latex');
  ylabel(['Principle axis ',num2str(2)],'interpreter', 'latex');
  zlabel(['Principle axis ',num2str(3)],'interpreter', 'latex');
end
end