
clear all; close all; clc;
mtP = [0.1;0.05];
lipW = [0.1;-0.15];
lipTh = 50;
[time,mt1,mt2,lip,RT] = lip_activity(mtP,lipW,lipTh);
w = 0.2;
ind = find(lip==1);
tTmp = time(ind);
tmp = lip(ind);
p1 = plot([tTmp;tTmp],[1.4+tmp+w;1.4+tmp-w],'k');
hold on
ind = find(mt2==1);
tTmp = time(ind);
tmp = mt2(ind);
p2 = plot([tTmp;tTmp],[1+tmp+w;1+tmp-w],'k');
ind = find(mt1==1);
tTmp = time(ind);
tmp = mt1(ind);
p3 = plot([tTmp;tTmp],[0.6+tmp+w;0.6+tmp-w],'k');
axis tight
legend([p1(1),p2(1),p3(1)],"LIP","MT Inhibitory","MT Excitatory")
title("Raster plot")
xlabel("time")
ylabel("Spikes")


function [time,mt1,mt2,lip,RT] = lip_activity(mtP,lipW,lipTh)
    dt = 0.001;
    t = 0;
    N = [0;0];
    rate = 0;
    mt1 = [];
    mt2 = [];
    lip = [];
    time = [];
    lipT = [];
    M = 100;
    while rate < lipTh
        time = [time,t];
        dN = rand(2,1) < mtP;
        mt1 = [mt1,dN(1)];
        mt2 = [mt2,dN(2)];
        N = N + dN;
        p_lip = sum(N .* lipW);
        lipEvent = rand()<p_lip;
        lip = [lip,lipEvent];
        if (lipEvent == 1)
            lipT = [lipT,t];
        end
        if (length(lipT) >= M)
            rate = M / (t-lipT(end-M+1));
        end
        t = t + dt;
    end
    RT = t;
end