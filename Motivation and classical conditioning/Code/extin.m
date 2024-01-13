clc
clear
trial=200;
Lrate=0.05;
reward=[];
reward=[reward ones(1,100)];
reward=[reward zeros(1,100)];
u =1;
w = zeros(1,200);
for i = 1:199
     V=u*w(1,i);
     d= reward(i)-V;
     w(1,i+1)=w(1,i)+Lrate*d*u;
end
t=1:trial;
scatter(t,w)
title("Extinction")
xlabel("trials",'Fontsize',16)
ylabel("w",'Fontsize',16)
xline(100,'--r',{'Start','train'})
xline(0,'--r',{'Pre-','train'})
