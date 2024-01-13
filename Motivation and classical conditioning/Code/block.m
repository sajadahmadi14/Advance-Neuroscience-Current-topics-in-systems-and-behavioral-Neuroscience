clc
clear
trial=200;
Lrate=0.05;

% u=ones(1,200);
% u2=[];
% u2=[u2 zeros(1,100)];
% u2=[u2 ones(1,100)];
% %u2=zeros(1,200);
u = [ones(1,trial);zeros(1,trial/2),ones(1,trial/2)];
reward=ones(1,200);
w = zeros(2,200);
%w2 = zeros(1,200);
    for i = 1:trial-1
        delta = reward(i) - u(:,i)' * w(:,i);
        w(:,i+1) = w(:,i) + Lrate*delta*u(:,i);
    end

t=1:trial;
scatter(t,w(1,:))
hold on
scatter(t,w(2,:))
title("Blocking")
xlabel("trials",'Fontsize',16)
ylabel("w",'Fontsize',16)
xline(100,'--r',{'Start','train'})
xline(0,'--r',{'Pre-','train'})
