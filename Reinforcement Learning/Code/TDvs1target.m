clear
close all
clc

lambda=0.95;
DF=0.01:0.25:1;
DF=[DF 1];
LR=0.01:0.25:1;
LR=[LR 1];
punishloc = [5,5];
rewardloc = [10,10];
rewardval=50;
punishval=-50;

vmap=zeros(15,15);
vmap(punishloc(1),punishloc(2))=punishval;
vmap(rewardloc(1),rewardloc(2))=rewardval;


%% main
for ii=1:size(LR,2)
    for jj=1:size(DF,2)
        lr=LR(ii);
        df=DF(jj);
        pmap = ones(15,15,4);
        for tt=1:100
            ss=1;
            clear agent
            if tt==100
                agent(ss,:)=[2,13];
            else
                agent(ss,:)=randi([1,15],[1,2]);
            end
            
            criteria=0;
            while criteria==0
                action=act(agent(ss,:),pmap);
                ACTION(ss) = action;
                agent(ss+1,:)=move(agent(ss,:),action);
                
                reward = vmap(agent(ss+1,1),agent(ss+1,2));
                ss=ss+1;
                if reward ~= 0
                    criteria=1;
                end
                
            end
            
            steps = size(agent,1);
            for zz=steps-1:-1:1
                
                
                d=reward + df*(mean(pmap(agent(steps,1),agent(steps,2),ACTION(zz))) - mean(pmap(agent(steps,1),agent(steps,2),ACTION(zz))));
                d=d*((lambda)^(steps-zz-1));
                pmap(agent(zz,1),agent(zz,2),ACTION(zz))=pmap(agent(zz,1),agent(zz,2),ACTION(zz))+d*lr;
                
            end
            Pmap{tt}=mean(pmap,3);
            stepnum(tt)=ss;
            AGENT{tt}=agent;
        end
        STEPnum(ii,jj,:)=stepnum;
        Apath{ii,jj}=AGENT(100);
        Pmap{ii,jj}=mean(pmap,3);
        clear pmap ss d reward action agent tt ACTION AGENT stepnum steps
    end
end

%% number of steps
figure(1)
imagesc(STEPnum(:,:,100))
colorbar
set(gca, 'YDir','reverse')
subtitle('(Number of Steps)','fontweight','bold','fontsize',15)
title('Learning Rate vs. Discount Factor','fontweight','bold','fontsize',15)
ylabel('Learning Rate','fontweight','bold')
xlabel('Discount Factor','fontweight','bold')
xticks([1 2 3 4 5])
xticklabels({'0.01','0.26','0.51','0.76','1'})
yticks([1 2 3 4 5])
yticklabels({'0.01','0.26','0.51','0.76','1'})

%% path
c1=0;
for ii=1:5
    for jj=1:5
        c1=c1+1;;
        figure(2)
        subplot(5,5,c1)
        imagesc(vmap)
        hold on
        agent=cell2mat(Apath{ii, jj});
        for kk=1:size(agent,1)-1
            p1=agent(kk,:);
            p2=agent(kk+1,:);
            dp=p2-p1;
            subplot(5,5,c1)
            quiver(p1(2),p1(1),dp(2),dp(1),0,'w')
        end
        clear agent
    end
end

figure(2)
sgtitle('          Learning Rate vs. Discount Factor','fontweight','bold','fontsize',15)
subplot(5,5,3); title('(Agents Path)','fontweight','bold','fontsize',15)

subplot(5,5,1); ylabel('LR=0.1','fontweight','bold')
subplot(5,5,6); ylabel('LR=0.26','fontweight','bold')
subplot(5,5,11); ylabel('LR=0.51','fontweight','bold')
subplot(5,5,16); ylabel('LR=0.76','fontweight','bold')
subplot(5,5,21); ylabel('LR=1','fontweight','bold')

subplot(5,5,21); xlabel('DF=0.1','fontweight','bold')
subplot(5,5,22); xlabel('DF=0.26','fontweight','bold')
subplot(5,5,23); xlabel('DF=0.51','fontweight','bold')
subplot(5,5,24); xlabel('DF=0.76','fontweight','bold')
subplot(5,5,25); xlabel('DF=1','fontweight','bold')

%% pmap
c1=0;
for ii=1:5
    for jj=1:5
        c1=c1+1;;
        figure(3)
        subplot(5,5,c1)
        hold on
        agent=imagesc((Pmap{ii, jj}));
set(gca, 'Ydir', 'reverse');
        clear agent
    end
end

figure(3)
sgtitle('          Learning Rate vs. Discount Factor','fontweight','bold','fontsize',15)
subplot(5,5,3); title('(State Action Value)','fontweight','bold','fontsize',15)

subplot(5,5,1); ylabel('LR=0.1','fontweight','bold')
subplot(5,5,6); ylabel('LR=0.26','fontweight','bold')
subplot(5,5,11); ylabel('LR=0.51','fontweight','bold')
subplot(5,5,16); ylabel('LR=0.76','fontweight','bold')
subplot(5,5,21); ylabel('LR=1','fontweight','bold')

subplot(5,5,21); xlabel('DF=0.1','fontweight','bold')
subplot(5,5,22); xlabel('DF=0.26','fontweight','bold')
subplot(5,5,23); xlabel('DF=0.51','fontweight','bold')
subplot(5,5,24); xlabel('DF=0.76','fontweight','bold')
subplot(5,5,25); xlabel('DF=1','fontweight','bold')










%% function

function nextloc=move(agent,action)
direction(1,:)=[0,-1];
direction(2,:)=[-1,0];
direction(3,:)=[0,1];
direction(4,:)=[1,0];

nextloc=agent+direction(action,:);

if nextloc(1) > 15
    nextloc(1) = 15;
end

if nextloc(1) < 1
    nextloc(1) = 1;
end

if nextloc(2) > 15
    nextloc(2) = 15;
end

if nextloc(2) < 1
    nextloc(2) = 1;
end
end

function action=act(agent,pmap)
action=wheel(pmap(agent(1),agent(2),:));
end

function action=wheel(probabilities)
probabilities = probabilities - min(probabilities)+1;
x = cumsum(probabilities);
x = x./x(4);
n = rand(1,1);

if n<x(1)
    action = 1;
elseif n<x(2)
    action = 2;
elseif n<x(3)
    action = 3;
else
    action = 4;
end
end
