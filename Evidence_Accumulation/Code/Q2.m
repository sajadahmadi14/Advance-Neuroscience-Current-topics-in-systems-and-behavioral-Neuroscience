clc
clear all
close all
T=1;
dt=0.1;
x=cell(1,20);
c=zeros(1,20);
for l=1:20
[c1 x1]=sm(1,1,0,dt,0,T);
c(1,l)=c1;
x{1, l}  =x1;
clear x1 c1
hold on
scatter(l,c(1,l),50,'filled');

xlabel('Choice')
ylabel('trials')
end
figure;
for p=1:20
subplot(4,5,p)
hold on
if c(1,p)==-1
plot([0: dt: T-dt], x{1,p},'r');
else
    plot([0: dt: T-dt], x{1,p});
end
title('Experiment',num2str(p))
xlabel('time')
ylabel('x')
end

B=[1,10,0,0.1,-1];
T=10;
dt=0.1;
figure;
for w=1:6
    if w==6
        for q=1:5
       subplot(2,3,w) 
       B=[1,10,0,0.1,-1];
B=B(q);
[c2 x2]=sm(B,1,0,dt,0,T);
plot([0: dt: T-dt], x2, 'linewidth', 1.5);
hold on
legend('B = 1', 'B = 10', 'B = 0', 'B = 0.1', 'B = -1');
xlabel('time')
title('Different Bias values')
ylabel('x')
        end
else
subplot(2,3,w)
B=[1,10,0,0.1,-1];
B=B(w);
[c2 x2]=sm(B,1,0,dt,0,T);
plot([0: dt: T-dt], x2, 'linewidth', 1.5);
hold on
legend(['B = ' num2str(B)]);

xlabel('time')
title('Different Bias values')
ylabel('x')
        end
end

function [choice, x] = sm(B, sigma, mean, dt, x0, T)
    x(1) = x0;
    t=T/dt;
    for i = 2:t
        x(i) = x(i-1) + B*dt + sigma * normrnd(mean, sqrt(dt)); 
    end
    choice = sign(x(:,T/dt));
end

 

 