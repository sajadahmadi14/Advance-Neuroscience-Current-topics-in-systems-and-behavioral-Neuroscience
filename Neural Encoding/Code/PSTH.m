
clc
clear all
load UnitsData.mat
n=385;
h =[];
h0=[];
h1=[];
h2=[];
h3=[];
h4=[];
h5=[];
h6=[];

index=Unit(n).Cnd(1).TrialIdx;
A=Unit(n).Trls(index) ;  
B=Unit(n).Trls;  

for k=1:length(B)
h(k,:)=hist(B{k,1},53);
end
h=mean(h);


for k=1:length(A)
h0(k,:)=hist(A{k,1},53);
end
h0=mean(h0);

index=Unit(n).Cnd(2).TrialIdx;
A=Unit(n).Trls(index) ;  
for k=1:length(A)
h1(k,:)=hist(A{k,1},53);
end
h1=mean(h1);
index=Unit(n).Cnd(3).TrialIdx;
A=Unit(n).Trls(index) ;  
for k=1:length(A)
h2(k,:)=hist(A{k,1},53);
end
h2=mean(h2);
index=Unit(n).Cnd(4).TrialIdx;
A=Unit(n).Trls(index) ;  
for k=1:length(A)
h3(k,:)=hist(A{k,1},53);
end
h3=mean(h3);
index=Unit(n).Cnd(5).TrialIdx;
A=Unit(n).Trls(index) ;  
for k=1:length(A)
h4(k,:)=hist(A{k,1},53);
end
h4=mean(h4);
index=Unit(n).Cnd(6).TrialIdx;
A=Unit(n).Trls(index) ;  
for k=1:length(A)
h5(k,:)=hist(A{k,1},53);
end
h5=h5*5000;
h5=mean(h5);
%subplot(2,2,4);
t=linspace(-1.2,2,length(h));
plot(t,h,':','Linewidth',2)
%p(1).LineWidth = 100;
xlim([-1.2 2])
ylim([0 max(h)+1])
hold on
plot(t,h0)
plot(t,h1)
plot(t,h2)
plot(t,h3)
plot(t,h4)
plot(t,h5)
plot([0 0],[0 max(h)+1])
xlabel('Time');
ylabel('spike');
title(['PSTH of unit ' ,num2str(n), ' confitions']);
legend('All trials','[3 -1]','[3 1]','[6 -1]',...
    '[6 1]','[9 -1]','[9 1]');

