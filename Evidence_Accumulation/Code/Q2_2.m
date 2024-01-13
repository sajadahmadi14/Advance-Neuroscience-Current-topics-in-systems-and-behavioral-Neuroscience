clear all; close all; clc;

act1 = @(x) sin(x+pi/4)^10;
act2 = @(y) sin(y-pi/4)^10;
dt = 0.001;
T = 0.5;
L = T/dt;
lipW = [0.1;-0.1];
mtTheta = [[ones(1,200)*0.2,ones(1,50)*0.8,ones(1,250)*0.2] ; ones(1,L)*1.85];
[time,mt1,mt2,lip1,lip2] = lip_activity2(mtTheta,act1,act2,lipW,T);
w = 0.2;
ind = find(lip1==1);
tTmp = time(ind);
tmp = lip1(ind);
p1 = plot([tTmp;tTmp],[2+tmp+w;2+tmp-w],'k');
ind = find(lip2==1);
tTmp = time(ind);
tmp = lip2(ind);
hold on
p2 = plot([tTmp;tTmp],[1.5+tmp+w;1.5+tmp-w],'k');
ind = find(mt1==1);
tTmp = time(ind);
tmp = mt1(ind);
p3 = plot([tTmp;tTmp],[1+tmp+w;1+tmp-w],'k');
ind = find(mt2==1);
tTmp = time(ind);
tmp = mt2(ind);
p4 = plot([tTmp;tTmp],[0.5+tmp+w;0.5+tmp-w],'k');
xlim([0,T])
axis tight
title("Raster plot")
xlabel("time")
ylabel("Spike")
legend([p1(1),p2(1),p3(1),p4(1)],"LIP1 = +MT1 - MT2","LIP2 = -MT1 + MT2","MT1","MT2")


function [time,mt1,mt2,lip1,lip2] = lip_activity2(mtTheta,act1,act2,lipW,T)
    dt = 0.001;
    t = 0;
    N = [0;0];
    mt1 = [];
    mt2 = [];
    lip1 = [];
    lip2 = [];
    time = [];
    cnt = 1;
    lipW1 = lipW;
    lipW2 = flip(lipW);
    for i = 1:round(T/dt)
        time = [time,t];
        
        theta = mtTheta(:,cnt);
        dN = rand(2,1) < [act1(theta(1));act2(theta(2))];
        mt1 = [mt1,dN(1)];
        mt2 = [mt2,dN(2)];
        N = N + dN;
        p_lip1 = sum(N .* lipW1);
        p_lip2 = sum(N .* lipW2);
        lipEvent1 = rand(1)<p_lip1;
        lipEvent2 = rand(1)<p_lip2;
        lip1 = [lip1,lipEvent1];
        lip2 = [lip2,lipEvent2];
        t = t + dt;
        cnt = cnt+1;
    end
end