clear
close all
clc

df=0.5;
lr=0.5;
punishloc = [5,5];
rewardloc = [10,10];
rewardval=5;
punishval=-5;

vmap=zeros(15,15);
vmap(punishloc(1),punishloc(2))=punishval;
vmap(rewardloc(1),rewardloc(2))=rewardval;

pmap = ones(15,15,4);

%% main
for tt=1:1000
    ss=1;
    clear agent
    if ismember(tt,[1 10 250 500 1000])
        agent(ss,:)=[1,13];
    else
    	agent(ss,:)=randi([1,15],[1,2]);
    end
    
    criteria=0;
    while criteria==0
        action=act(agent(ss,:),pmap);
        agent(ss+1,:)=move(agent(ss,:),action);
        
        reward = vmap(agent(ss+1,1),agent(ss+1,2));
        d=reward + df*(mean(pmap(agent(ss+1,1),agent(ss+1,2),:)) - mean(pmap(agent(ss,1),agent(ss,2),:)));
        pmap(agent(ss,1),agent(ss,2),action) = pmap(agent(ss,1),agent(ss,2),action) + d* lr;
        ss=ss+1;
        if reward ~= 0
            criteria=1;
        end
        
        
%-----------------------movement of agent each step-----------------------%        
        
% figure(1)
% imagesc(vmap)
% hold on
% scatter(agent(ss,2),agent(ss,1),100,'k','filled','linewidth',5)
% scatter(punishloc(2),punishloc(1),100,'r','linewidth',5)
% scatter(rewardloc(2),rewardloc(1),100,'g','linewidth',5)
% 
% title(num2str(tt))
% subtitle(num2str(ss))
% pause(0.001)
    end
    
    Pmap{tt}=mean(pmap,3);
    stepnum(tt)=ss;
    AGENT{tt}=agent;
end
%% agents path in each chosen trial
trial=500;

figure(98)
sgtitle(['Trial = ',num2str(trial)])
subplot(2,4,[1,2,5,6])
imagesc(vmap)
hold on

agent=cell2mat(AGENT(trial));
for ii=1:size(agent,1)-1
    p1=agent(ii,:);
    p2=agent(ii+1,:);
    dp=p2-p1;
    subplot(2,4,[1,2,5,6])
    quiver(p1(2),p1(1),dp(2),dp(1),0,'w')
%     pause(0.001)
end
scatter(agent(1,2),agent(1,1),100,'filled','k','linewidth',5)
[px,py]=gradient(cell2mat(Pmap(trial)));
subplot(2,4,[3,4,7,8])
imagesc(cell2mat(Pmap(trial)))
hold on
quiver(px,py,'w')
scatter(punishloc(2),punishloc(1),100,'r','linewidth',5)
scatter(rewardloc(2),rewardloc(1),100,'g','linewidth',5)
set(gca, 'YDir','reverse')


subplot(2,4,[1,2,5,6])
title('Agents Path')
subplot(2,4,[3,4,7,8])
title('State Action Value')


%% demo
trial=1000;
agent=cell2mat(AGENT(trial));


writerObj = VideoWriter('Demo.avi');
writerObj.FrameRate = 10;
open(writerObj);
for ss=1:size(agent,1)-1
figure(1)
imagesc(vmap)
hold on
scatter(agent(ss,2),agent(ss,1),100,'k','filled','linewidth',5)
scatter(punishloc(2),punishloc(1),100,'r','linewidth',5)
scatter(rewardloc(2),rewardloc(1),100,'g','linewidth',5)
scatter(agent(1,2),agent(1,1),100,'filled','w','linewidth',5)


title(['Trial = ',num2str(trial)])
subtitle(['# Steps = ',num2str(ss)])
% pause(0.01)



    F = getframe(gcf) ;
    drawnow
    writeVideo(writerObj, F);
end
close(writerObj);





























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
