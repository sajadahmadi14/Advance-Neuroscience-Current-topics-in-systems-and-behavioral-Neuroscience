clc
clear
trial=400;
Lrate=0.05;

%u = [ones(1,trial) ; rand(1,trial) < 0.2];
u=zeros(2,trial);
u(1,:)=ones(1,trial);
a=randi([0:1], [1,400]);
u(2,:)=a;
reward = 1-u(2,:);
    w = zeros(2,trial);
    w(:,1) = 0;
    for i = 1:trial-1
        delta = reward(i) - u(:,i)' * w(:,i);
        w(:,i+1) = w(:,i) + Lrate*delta*u(:,i);
    end

t=1:trial;
scatter(t,w)
title("Inhinitory")
xlabel("trials",'Fontsize',16)
ylabel("w",'Fontsize',16)
xline(100,'--r',{'Start','train'})
xline(0,'--r',{'Pre-','train'})

