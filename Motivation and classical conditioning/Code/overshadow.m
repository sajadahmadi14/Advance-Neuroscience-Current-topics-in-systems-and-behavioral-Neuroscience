clc
clear
nTrials=200;
etha=0.05;

u = zeros(2,nTrials);
u(1,:)=rand(1,nTrials)<0.6;
u(2,:)=rand(1,nTrials)<0.4;
r = ones(1,nTrials);
w=zeros(2,nTrials);
   for i = 2:nTrials-1
        delta = r(i) - u(:,i)' * w(:,i);
        w(:,i+1) = w(:,i) + etha*delta*u(:,i);
   end
scatter(1:nTrials,w(1,:))
hold on
scatter(1:nTrials,w(2,:))
ylim([0 1])
title("Overshadow")
xlabel("trials",'Fontsize',16)
ylabel("w",'Fontsize',16)


w1=mean(w(1,:));
w2=mean(w(2,:));
yline(w1,'-','mean','LineWidth',3);
yline(w2,'-','mean','LineWidth',3);
 

