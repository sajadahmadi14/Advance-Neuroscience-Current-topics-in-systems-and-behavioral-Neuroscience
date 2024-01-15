%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_6



%% part 1&2
close all;clc;clear;
r=zeros(15,15); 
x1=9 ; y1=10;
x2=8 ; y2=7;
r(x1,y1)=20;
r(x2,y2)=-20; 
direction=[1,2,3,4]; % 1:left , 2:Right , 3:bottom , 4:up
Q=Qfunc();
Trials = 450;
eta = 0.6;
gamma = 0.4;
for nTrial = 1:Trials
rng shuffle
i=randi(15); 
j=randi(15); 
if i==x1 
   yTemp=[1:y1-1,y1+1:15];
   c=randi(14);
   j=yTemp(c);
end
if i==x2 
   yTemp=[1:y2-1,y2+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[1:x1-1,x1+1:15];
   c=randi(14);
   i=xTemp(c);
end
if j==y2 
   xTemp=[1:x2-1,x2+1:15];
   c=randi(14);
   i=xTemp(c);
end   
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end
        iLast = i;
        jLast = j;
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        Step = Step + 1;
    end
    PosStruct{nTrial} = posVec;
    mStep(nTrial)=abs(posVec(1,1)-x1)+abs(posVec(1,2)+1-y1);
    for p=1:size(posVec,2)
    path(nTrial).pmat=zeros(15,15);
    end
    for p=1:size(posVec,2)
        path(nTrial).pmat(posVec(1,p),posVec(2,p))=50;
    end
end

for stepc=1:Trials
   stepcounter(stepc)=length(Qstruct{1,stepc}.Step);
end

for mm=1:Trials-7
   n(mm)= length(find(stepcounter(mm:mm+4)<mStep(nTrial) + 6));
    if n(mm) > 4
          Tnew = mm;
          break;
    end
end

figure;
Trial=[10,15,20,Tnew,Tnew+1,Tnew+2,Tnew+3,Tnew+4];
for i = 1:length(Trial)
    f=figure;
    colormap(white)
    posVec = PosStruct{Trial(i)};
    maze=zeros(15,15);
    imagesc(maze,'AlphaData',0.2)
    colormap(winter)
    h = gca;  
    set(h, 'YDir', 'normal');
    for x=1:15
        for y=1:15
    if path(Trial(i)).pmat(x,y)==50
       text(x-0.4,y-0.25,'\bullet','Color','blue','FontSize',35)
    end
        end
    end
    text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
    text(x2-0.25,y2,'X','Color','red','FontSize',15)
    text(posVec(1,1)-0.3,posVec(2,1)-0.15,'\bullet','Color','yellow','FontSize',25)
    xlim([1,15])
    ylim([1,15])
    title(['Trial ',num2str(Trial(i))],'interpreter', 'latex')  
     fig_name = strcat('part10_',num2str(i));
     print(f,fig_name,'-depsc'); 
end

figure;
colormap(jet)
pcolor(1:1:15,1:1:15,log10(max(Q,[],3))')
h = gca;  
colorbar
text(x1-0.3,y1+0.15,'\bullet','Color','green','FontSize',25)
text(x2-0.25,y2,'X','Color','red','FontSize',15)
print -depsc part1_a.eps

figure;
colormap(jet)
Q(x1,y1,:)=20;
Q(x1,y1,:)=1;
contourf(1:1:15,1:1:15,log10(max(Q,[],3))')
hold on
text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
text(x2-0.25,y2,'X','Color','red','FontSize',15)
colorbar
axis square
print -depsc part1_b.eps

figure;
[fx,fy] = gradient(max(Q,[],3)');
quiver(1:15,1:15,fx,fy,'k','AutoScaleFactor',0.6);
hold on
text(x1,y1,'\bullet','Color','green','FontSize',25)
text(x2,y2,'X','Color','red','FontSize',15)
xlim([1 15])
ylim([1 15])
print -depsc part1_c.eps

figure;
set(gcf,'units','points','position',[0,0,600,350])
v = VideoWriter('maze.avi');
open(v)
maze=zeros(15,15);
imagesc(maze,'AlphaData',0.2)
colormap(winter)
h = gca;  
set(h, 'YDir', 'normal');
for nTrial = [10,15,20,Tnew,Tnew+1,Tnew+2,Tnew+3,Tnew+4]
    Path = PosStruct{nTrial};
    x=[];
    y=[];
    x = Path(1,:);
    x=[x x1];
    y = Path(2,:);
    y=[y y1];
    for i = 1:length(Path)-1
        step=i;
        scatter(x1-0.35,y1+0.2,100,'g','filled')
        hold on
        scatter(x2-0.35,y2+0.2,150,'rX')
        hold on
        scatter(x(1)-0.35,y(1)+0.2,100,'y','filled')
        hold on
        scatter(x(i)-0.35,y(i)+0.2,100,'blue','filled')
        hold off
        xlim([0 16])
        ylim([0 16])
        title(sprintf("Trial %d , Step= %d",nTrial,step),'interpreter','latex')
        frame = getframe(gcf);
        for i = 1:5
            writeVideo(v,frame);
        end
    end
end
close(v)

figure;
set(gcf,'units','points','position',[0,0,600,350])
v = VideoWriter('value1.avi');
open(v)
for  nTrial = [1:100,150,200,250,300,350,400,450]
        Q=Qstruct{1,nTrial}.Step{1,length(Qstruct{1,nTrial}.Step)};
        colormap(jet)
        Q(x1,y1,:)=20;
        Q(x1,y1,:)=1;
        contourf(1:1:15,1:1:15,log10(max(Q,[],3))')
        hold on
        text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
        text(x2-0.25,y2,'X','Color','red','FontSize',15)
        hold off
        colorbar
        xlim([0 16])
        ylim([0 16])
        title(sprintf("Trial %d " ,nTrial),'interpreter','latex')
        frame = getframe(gcf);
        for i = 1:5
            writeVideo(v,frame);
        end
    end
close(v)

figure;
set(gcf,'units','points','position',[0,0,600,350])
v = VideoWriter('value2.avi');
open(v)
for nTrial = [1:100,150,200,250,300,350,400,450]
        Q=Qstruct{1,nTrial}.Step{1,length(Qstruct{1,nTrial}.Step)};
        colormap(jet)
        pcolor(1:1:15,1:1:15,log10(max(Q,[],3))')
        h = gca;  
        colorbar
        text(x1-0.3,y1+0.15,'\bullet','Color','green','FontSize',25)
        text(x2-0.25,y2,'X','Color','red','FontSize',15)
        hold off
        xlim([0 16])
        ylim([0 16])
        title(sprintf("Trial %d " ,nTrial),'interpreter','latex')
        frame = getframe(gcf);
        for i = 1:5
            writeVideo(v,frame);
        end
    end
close(v)

figure;
set(gcf,'units','points','position',[0,0,600,350])
v = VideoWriter('gradients.avi');
open(v)
for nTrial = [1:100,150,200,250,300,350,400,450]
        Q=Qstruct{1,nTrial}.Step{1,length(Qstruct{1,nTrial}.Step)};
        [fx,fy] = gradient(max(Q,[],3)');
        quiver(1:15,1:15,fx,fy,'k','AutoScaleFactor',0.6);
        hold on
        text(x1,y1,'\bullet','Color','green','FontSize',25)
        text(x2,y2,'X','Color','red','FontSize',15)
        hold off
        xlim([0 16])
        ylim([0 16])
        title(sprintf("Trial %d " ,nTrial),'interpreter','latex')
        frame = getframe(gcf);
        for i = 1:5
            writeVideo(v,frame);
        end
    end
close(v)
%% part 3
close all;clc;clear;
gammaVec=[0.1:0.1:1];
etaVec=[0.1:0.1:1];
for ev=1:length(etaVec)
    eta=etaVec(ev);
    for gv=1:length(gammaVec)
        gamma=gammaVec(gv)
for iter=1:10
r=zeros(15,15); 
x1=3 ; y1=10;
x2=8 ; y2=7;
r(x1,y1)=10; 
r(x2,y2)=-10; 
direction=[1,2,3,4]; 
Q=Qfunc();
Trials = 450;
for nTrial = 1:Trials
rng shuffle
i=randi(15); 
j=randi(15); 
if i==x1 
   yTemp=[1:y1-1,y1+1:15];
   c=randi(14);
   j=yTemp(c);
end
if i==x2 
   yTemp=[1:y2-1,y2+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[1:x1-1,x1+1:15];
   c=randi(14);
   i=xTemp(c);
end
if j==y2 
   xTemp=[1:x2-1,x2+1:15];
   c=randi(14);
   i=xTemp(c);
end   
    i=8;j=10;
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end     
        iLast = i;
        jLast = j;
        
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        Step = Step + 1;
    end
    PosStruct{nTrial} = posVec;
    mStep(nTrial)=abs(posVec(1,1)-x1)+abs(posVec(1,2)+1-y1);
    for p=1:size(posVec,2)
    path(nTrial).pmat=zeros(15,15);
    end
    for p=1:size(posVec,2)
        path(nTrial).pmat(posVec(1,p),posVec(2,p))=50;
    end
end
for stepc=1:Trials
   stepcounter(stepc)=length(Qstruct{1,stepc}.Step);
end
for mm=1:Trials-7
   n(mm)= length(find(stepcounter(mm:mm+4)<mStep(nTrial) + 6));
    if n(mm) > 4
          Tnew = mm;
          break;
    end
end
save('Tnew.mat','Tnew');
clear stepcounter
clear Qstruct
end
    end
end

load('Tnew.mat')
f=median(Tnew,1)
subplot(2,1,1)
plot(0.1:0.1:1,f(:,:,4),'b')
xlabel('$\epsilon(Learning Rate)$','interpreter','latex')
ylabel('Minimum number of Trials','interpreter','latex')
title('$\gamma=0.4$','interpreter','latex')
subplot(2,1,2)
ff=[];
for z=1:10
ff=[ff f(:,2,z)]
end
plot(0.1:0.1:1,ff,'b')
xlabel('$\gamma$','interpreter','latex')
ylabel('Minimum number of Trials','interpreter','latex')
title('$\epsilon=0.2$','interpreter','latex')
print -depsc part3.eps

%% part4_a1
close all;clc;clear;
gammaVec=[0.5];
etaVec=[0.1:0.1:1];
t1=zeros(10,450,10);
t2=zeros(10,450,10);
for ev=1:length(etaVec)
    eta=etaVec(ev);
    for gv=1:length(gammaVec)
        gamma=gammaVec(gv)
for iter=1:10
r=zeros(15,15); 
x1=3 ; y1=8;
x2=13 ; y2=8;
r(x1,y1)=30; 
r(x2,y2)=5; 
direction=[1,2,3,4];
Q=Qfunc();
Trials = 450;
for nTrial = 1:Trials
rng shuffle
i=randi(15); 
j=randi(15);
if i==x1 
   yTemp=[1:y1-1,y1+1:15];
   c=randi(14);
   j=yTemp(c);
end
if i==x2 
   yTemp=[1:y2-1,y2+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[1:x1-1,x1+1:15];
   c=randi(14);
   i=xTemp(c);
end
if j==y2 
   xTemp=[1:x2-1,x2+1:15];
   c=randi(14);
   i=xTemp(c);
end   
    i=8;
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end              
        iLast = i;
        jLast = j;
        
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        if (i==x1 & j==y1)
            t1(iter,nTrial,ev)=t1(iter,nTrial,ev)+1;
        end
        if (i==x2 & j==y2)
            t2(iter,nTrial,ev)=t2(iter,nTrial,ev)+1; 
        end
        save('t1_gamma.mat','t1');
        Step = Step + 1;
    end
end
clear stepcounter
clear Qstruct
end
    end
end
%% part4_a2
close all;clc;clear;
gammaVec=[0.1:0.1:1];
etaVec=[0.4];
t1=zeros(10,450,10);
t2=zeros(10,450,10);
for ev=1:length(etaVec)
    eta=etaVec(ev);
    for gv=1:length(gammaVec)
        gamma=gammaVec(gv)
for iter=1:1
r=zeros(15,15); 
x1=3 ; y1=8;
x2=13 ; y2=8;
r(x1,y1)=30; 
r(x2,y2)=5; 
direction=[1,2,3,4];
Q=Qfunc();
Trials = 450;
for nTrial = 1:Trials
rng shuffle
i=randi(15); 
j=randi(15);
if i==x1 
   yTemp=[1:y1-1,y1+1:15];
   c=randi(14);
   j=yTemp(c);
end
if i==x2 
   yTemp=[1:y2-1,y2+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[1:x1-1,x1+1:15];
   c=randi(14);
   i=xTemp(c);
end
if j==y2 
   xTemp=[1:x2-1,x2+1:15];
   c=randi(14);
   i=xTemp(c);
end   
    i=8;
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end              
        iLast = i;
        jLast = j;
        
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        if (i==x1 & j==y1)
            t1(iter,nTrial,gv)=t1(iter,nTrial,gv)+1;
        end
        if (i==x2 & j==y2)
            t2(iter,nTrial,gv)=t2(iter,nTrial,gv)+1; 
        end
        save('t1_eta.mat','t1');
        Step = Step + 1;
    end
end
clear stepcounter
clear Qstruct
end
    end
end
%% part4_b1
close all;clc;clear;
gammaVec=[0.5];
etaVec=[0.1:0.1:1];
t1=zeros(10,450,10);
t2=zeros(10,450,10);
for ev=1:length(etaVec)
    eta=etaVec(ev);
    for gv=1:length(gammaVec)
        gamma=gammaVec(gv)
for iter=1:10
r=zeros(15,15); 
x1=3 ; y1=8;
x2=13 ; y2=8;
r(x1,y1)=30; 
r(x2,y2)=5; 
direction=[1,2,3,4];
Q=Qfunc();
Trials = 450;
for nTrial = 1:Trials
rng shuffle
i=randi(7); 
j=randi(15); 
if i==x1 
   yTemp=[1:y1-1,y1+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[1:x1-1,x1+1:7];
   c=randi(6);
   i=xTemp(c);
end
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end              
        iLast = i;
        jLast = j;
        
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        if (i==x1 & j==y1)
            t1(iter,nTrial,ev)=t1(iter,nTrial,ev)+1;
        end
        if (i==x2 & j==y2)
            t2(iter,nTrial,ev)=t2(iter,nTrial,ev)+1; 
        end
        save('t1_gamma_closer.mat','t1');
        Step = Step + 1;
    end
end
clear stepcounter
clear Qstruct
end
    end
end
%% part4_b2
close all;clc;clear;
gammaVec=[0.1:0.1:1];
etaVec=[0.4];
t1=zeros(10,450,10);
t2=zeros(10,450,10);
for ev=1:length(etaVec)
    eta=etaVec(ev);
    for gv=1:length(gammaVec)
        gamma=gammaVec(gv)
for iter=1:10
r=zeros(15,15); 
x1=3 ; y1=8;
x2=13 ; y2=8;
r(x1,y1)=30; 
r(x2,y2)=5; 
direction=[1,2,3,4];
Q=Qfunc();
Trials = 450;
for nTrial = 1:Trials
rng shuffle
i=randi(7); 
j=randi(15); 
if i==x1 
   yTemp=[1:y1-1,y1+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[1:x1-1,x1+1:7];
   c=randi(6);
   i=xTemp(c);
end
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end              
        iLast = i;
        jLast = j;
        
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        if (i==x1 & j==y1)
            t1(iter,nTrial,gv)=t1(iter,nTrial,gv)+1;
        end
        if (i==x2 & j==y2)
            t2(iter,nTrial,gv)=t2(iter,nTrial,gv)+1; 
        end
        save('t1_eta_closer.mat','t1');
        Step = Step + 1;
    end
end
clear stepcounter
clear Qstruct
end
    end
end
%% part4_c1
close all;clc;clear;
gammaVec=[0.5];
etaVec=[0.1:0.1:1];
t1=zeros(10,450,10);
t2=zeros(10,450,10);
for ev=1:length(etaVec)
    eta=etaVec(ev);
    for gv=1:length(gammaVec)
        gamma=gammaVec(gv)
for iter=1:10
r=zeros(15,15); 
x1=3 ; y1=8;
x2=13 ; y2=8;
r(x1,y1)=30; 
r(x2,y2)=5; 
direction=[1,2,3,4];
Q=Qfunc();
Trials = 450;
for nTrial = 1:Trials
rng shuffle
i=randi(7); 
i=i+8;
j=randi(15); 
if i==x2 
   yTemp=[1:y2-1,y2+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[9:x2-1,x2+1:15];
   c=randi(6);
   i=xTemp(c);
end
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end              
        iLast = i;
        jLast = j;
        
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        if (i==x1 & j==y1)
            t1(iter,nTrial,ev)=t1(iter,nTrial,ev)+1;
        end
        if (i==x2 & j==y2)
            t2(iter,nTrial,ev)=t2(iter,nTrial,ev)+1; 
        end
        save('t1_gamma_notcloser.mat','t1');
        Step = Step + 1;
    end
end
clear stepcounter
clear Qstruct
end
    end
end
%% part4_c2
close all;clc;clear;
gammaVec=[0.1:0.1:1];
etaVec=[0.4];
t1=zeros(10,450,10);
t2=zeros(10,450,10);
for ev=1:length(etaVec)
    eta=etaVec(ev);
    for gv=1:length(gammaVec)
        gamma=gammaVec(gv)
for iter=1:10
r=zeros(15,15); 
x1=3 ; y1=8;
x2=13 ; y2=8;
r(x1,y1)=30; 
r(x2,y2)=5; 
direction=[1,2,3,4];
Q=Qfunc();
Trials = 450;
for nTrial = 1:Trials
rng shuffle
i=randi(7);
i=i+8;
j=randi(15);
if i==x2 
   yTemp=[1:y2-1,y2+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[9:x2-1,x2+1:15];
   c=randi(6);
   i=xTemp(c);
end
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end              
        iLast = i;
        jLast = j;
        
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta * (r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        Qstruct{nTrial}.Step{Step} = Q;
        
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        if (i==x1 & j==y1)
            t1(iter,nTrial,gv)=t1(iter,nTrial,gv)+1;
        end
        if (i==x2 & j==y2)
            t2(iter,nTrial,gv)=t2(iter,nTrial,gv)+1; 
        end
        save('t1_eta_notcloser.mat','t1');
        Step = Step + 1;
    end
end
clear stepcounter
clear Qstruct
end
    end
end

%% part4
close all;clc;clear;
load('t1_gamma.mat');
t1=sum(t1,2);
t1=mean(t1,1);
t1_gamma_equal=reshape(t1,1,10)/450;
load('t1_eta.mat');
t1=sum(t1,2);
t1=mean(t1,1);
t1_eta_equal=reshape(t1,1,10)/450;
load('t1_gamma_closer.mat');
t1=sum(t1,2);
t1=mean(t1,1);
t1_gamma_closer=reshape(t1,1,10)/450;
load('t1_eta_closer.mat');
t1=sum(t1,2);
t1=mean(t1,1);
t1_eta_closer=reshape(t1,1,10)/450;
load('t1_gamma_notcloser.mat');
t1=sum(t1,2);
t1=mean(t1,1);
t1_gamma_notcloser=reshape(t1,1,10)/450;
load('t1_eta_notcloser.mat');
t1=sum(t1,2);
t1=mean(t1,1);
t1_eta_notcloser=reshape(t1,1,10)/450;


x1=3 ; y1=8;
x2=13 ; y2=8;
figure;
set(gcf,'units','points','position',[0,0,1600,650])
    subplot(2,3,2)
    maze=zeros(15,15);
    text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
    text(x1-0.8,y1+0.3,'reward=30','Color','k','FontSize',8,'interpreter','latex')
    text(x2-0.3,y2-0.15,'\bullet','Color','blue','FontSize',25)
    text(x2-0.8,y2+0.3,'reward=5','Color','k','FontSize',8,'interpreter','latex')
    xlim([1,15])
    ylim([1,15])
    title('Starting point at the same distance from both targets','interpreter','latex')
    patch([7.5 8.5 8.5 7.5], [0 0 15 15], 'g','facealpha',0.3,'LineStyle','none')
    
subplot(2,3,4)
plot(0.1:0.1:1,t1_eta_equal)
xlim([0.1,1])
ylim([0,1])
xlabel('$\gamma$','interpreter','latex')
ylabel('$Ratio$','interpreter','latex')
title('$\epsilon$=0.4','interpreter','latex')
subplot(2,3,6)
plot(0.1:0.1:1,t1_gamma_equal)
xlim([0.1,1])
ylim([0,1])
xlabel('$\epsilon$','interpreter','latex')
ylabel('$Ratio$','interpreter','latex')
title('$\gamma$=0.5','interpreter','latex')
print -depsc part4_1.eps

figure;
set(gcf,'units','points','position',[0,0,1600,650])
    subplot(2,3,2)
    maze=zeros(15,15);
    text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
    text(x1-0.8,y1+0.3,'reward=30','Color','k','FontSize',8,'interpreter','latex')
    text(x2-0.3,y2-0.15,'\bullet','Color','blue','FontSize',25)
    text(x2-0.8,y2+0.3,'reward=5','Color','k','FontSize',8,'interpreter','latex')
    xlim([1,15])
    ylim([1,15])
    title('Starting point close to the green target','interpreter','latex')
    patch([0 7.5 7.5 0], [0 0 15 15], 'g','facealpha',0.3,'LineStyle','none')

    subplot(2,3,4)
plot(0.1:0.1:1,t1_eta_closer)
xlim([0.1,1])
ylim([0.5,1])
xlabel('$\gamma$','interpreter','latex')
ylabel('$Ratio$','interpreter','latex')
title('$\epsilon$=0.4','interpreter','latex')
subplot(2,3,6)
plot(0.1:0.1:1,t1_gamma_closer)
xlim([0.1,1])
ylim([0.5,1])
xlabel('$\epsilon$','interpreter','latex')
ylabel('$Ratio$','interpreter','latex')
title('$\gamma$=0.5','interpreter','latex')
print -depsc part4_2.eps

figure;
set(gcf,'units','points','position',[0,0,1600,650])
    subplot(2,3,2)
    maze=zeros(15,15);
    text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
    text(x1-0.8,y1+0.3,'reward=30','Color','k','FontSize',8,'interpreter','latex')
    text(x2-0.3,y2-0.15,'\bullet','Color','blue','FontSize',25)
    text(x2-0.8,y2+0.3,'reward=5','Color','k','FontSize',8,'interpreter','latex')
    xlim([1,15])
    ylim([1,15])
    title('Starting point close to the blue target','interpreter','latex')
    patch([8.5 15 15 8.5], [0 0 15 15], 'g','facealpha',0.3,'LineStyle','none')
    
subplot(2,3,4)
plot(0.1:0.1:1,t1_eta_notcloser)
xlim([0.1,1])
ylim([0,1])
xlabel('$\gamma$','interpreter','latex')
ylabel('$Ratio$','interpreter','latex')
title('$\epsilon$=0.4','interpreter','latex')
subplot(2,3,6)
plot(0.1:0.1:1,t1_gamma_notcloser)
xlim([0.1,1])
ylim([0,1])
xlabel('$\epsilon$','interpreter','latex')
ylabel('$Ratio$','interpreter','latex')
title('$\gamma$=0.5','interpreter','latex')
print -depsc part4_3.eps

%% TD lambda
close all;clc;clear;
r=zeros(15,15); 
x1=9 ; y1=10;
x2=8 ; y2=7;
r(x1,y1)=20; 
r(x2,y2)=-20; 
direction=[1,2,3,4]; 
Q=Qfunc();
Trials = 450;
eta = 0.6;
gamma = 0.4;
for nTrial = 1:Trials
rng shuffle
i=randi(15); 
j=randi(15); 
if i==x1 
   yTemp=[1:y1-1,y1+1:15];
   c=randi(14);
   j=yTemp(c);
end
if i==x2 
   yTemp=[1:y2-1,y2+1:15];
   c=randi(14);
   j=yTemp(c);
end
if j==y1 
   xTemp=[1:x1-1,x1+1:15];
   c=randi(14);
   i=xTemp(c);
end
if j==y2 
   xTemp=[1:x2-1,x2+1:15];
   c=randi(14);
   i=xTemp(c);
end   
    is = i;
    js = j;
    flag = 0;
    posVec = [i;j];
    iLast = 0;
    jLast = 0;
    iLast2 = 0;
    jLast2 = 0;
    iLast3 = 0;
    jLast3 = 0;
    
    Step=1;
    while(flag ~= 1)
         act = find(Q(i,j,:) == max(Q(i,j,find(~isnan(Q(i,j,:))))));
         if (length(act)>1)
             rng;
             act = act(randi(length(act)));
         end
         
        iLast3 = iLast2;
        jLast3 = jLast2;
        iLast2 = iLast;
        jLast2 = jLast;
        iLast = i;
        jLast = j;
if act==1
i=i+1;
j=j;
end
if act==2
   i=i-1;
   j=j;
end
if act==3
   i=i;
   j=j+1;
end
if act==4
   i=i;
   j=j-1;
end
        posVec = [posVec,[i;j]];
        Q(iLast,jLast,act) = Q(iLast,jLast,act) + eta *0.7*(r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast,jLast,act));
        if (Step >= 2 && iLast2~=iLast && jLast2~=jLast )
            Q(iLast2,jLast2,act) = Q(iLast2,jLast2,act) + eta*0.7^2*(r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast2,jLast2,act));
        end
        if (Step >= 3 && iLast3~=iLast2 && jLast3~=jLast2 )
            Q(iLast3,jLast3,act) = Q(iLast3,jLast3,act) + eta*0.7^3*(r(i,j) + gamma*max(Q(i,j,:)) - Q(iLast3,jLast3,act));
        end
        Qstruct{nTrial}.Step{Step} = Q;
        if (i==x1 & j==y1) || (i==x2 & j==y2)
            flag=true;
        end
        Step = Step + 1;
    end
    PosStruct{nTrial} = posVec;
    
    mStep(nTrial)=abs(posVec(1,1)-x1)+abs(posVec(1,2)+1-y1);
    for p=1:size(posVec,2)
    path(nTrial).pmat=zeros(15,15);
    end
    for p=1:size(posVec,2)
        path(nTrial).pmat(posVec(1,p),posVec(2,p))=50;
    end
end
for stepc=1:Trials
   stepcounter(stepc)=length(Qstruct{1,stepc}.Step);
end

for mm=1:Trials-7
   n(mm)= length(find(stepcounter(mm:mm+4)<mStep(nTrial) + 6));
    if n(mm) > 4
          Tnew = mm;
          break;
    end
end

figure;
Trial=[10,15,20,Tnew,Tnew+1,Tnew+2,Tnew+3,Tnew+4];
for i = 1:length(Trial)
    f=figure;
    colormap(white)
    posVec = PosStruct{Trial(i)};
    maze=zeros(15,15);
    imagesc(maze,'AlphaData',0.2)
    colormap(winter)
    h = gca;  
    set(h, 'YDir', 'normal');
    for x=1:15
        for y=1:15
    if path(Trial(i)).pmat(x,y)==50
       text(x-0.4,y-0.25,'\bullet','Color','blue','FontSize',35)
    end
        end
    end
    text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
    text(x2-0.25,y2,'X','Color','red','FontSize',15)
    text(posVec(1,1)-0.3,posVec(2,1)-0.15,'\bullet','Color','yellow','FontSize',25)
    xlim([1,15])
    ylim([1,15])
    title(['Trial ',num2str(Trial(i))],'interpreter', 'latex')  
    % fig_name = strcat('part50_',num2str(i));
    % print(f,fig_name,'-depsc'); 
end
figure;
colormap(jet)
pcolor(1:1:15,1:1:15,log10(max(Q,[],3))')
h = gca;  
colorbar
text(x1-0.3,y1+0.15,'\bullet','Color','green','FontSize',25)
text(x2-0.25,y2,'X','Color','red','FontSize',15)
%print -depsc part5_a.eps

figure;
colormap(jet)
Q(x1,y1,:)=20;
Q(x1,y1,:)=1;
contourf(1:1:15,1:1:15,log10(max(Q,[],3))')
hold on
text(x1-0.3,y1-0.15,'\bullet','Color','green','FontSize',25)
text(x2-0.25,y2,'X','Color','red','FontSize',15)
colorbar
axis square
%print -depsc part5_b.eps

figure;
[fx,fy] = gradient(max(Q,[],3)');
quiver(1:15,1:15,fx,fy,'k','AutoScaleFactor',0.6);
hold on
text(x1,y1,'\bullet','Color','green','FontSize',25)
text(x2,y2,'X','Color','red','FontSize',15)
xlim([1 15])
ylim([1 15])
%print -depsc part5_c.eps

function Q=Qfunc()
Q=zeros(15,15,4);
for i=1:15
    for j=1:15
        if (i==1 && j==1)
                Q(i,j,2) = nan;
                Q(i,j,4) = nan;
        end
        if (i==1 && j==15)
                Q(i,j,2) = nan;
                Q(i,j,3) = nan;
        end
        if (i==15 && j==1)
                Q(i,j,1) = nan;
                Q(i,j,4) = nan;
        end
        if (i==15 && j==15)
                Q(i,j,1) = nan;
                Q(i,j,3) = nan;
        end
        if (i>=2 && i<=14 && j == 1)
                Q(i,j,4) = nan;
        end
        if (i>=2 && i<=14 && j == 15)
                Q(i,j,3) = nan;
        end
        if (j>=2 && j<=14 && i == 1)
                Q(i,j,2) = nan;
        end
        if (j>=2 && j<=14 && i == 15)
                Q(i,j,1) = nan;
        end
    end
end
end


