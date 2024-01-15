%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_3


%% Part 1
clear;clc;close all;

    filenames{1} = 'S_monkey1.mat';
    filenames{2} = 'S_monkey2.mat';
    filenames{3} = 'S_monkey3.mat';
    
    keepNeuronsName{1} = 'keepNeurons1.mat';
    keepNeuronsName{2} = 'keepNeurons2.mat';
    keepNeuronsName{3} = 'keepNeurons3.mat';
 
    monkeys = {'monkey1', 'monkey2', 'monkey3'};
    
    % Plot PSTH for Most Active Neurons in different Monkeys.
    for imonkey = 1:length(monkeys)
        load(filenames{imonkey});
        load(keepNeuronsName{imonkey});
        NeuronID=find(keepNeurons==1);
        num_gratings=length(S)
        f=figure;
        set(gcf,'units','points','position',[0,0,1200,600])
        for igrat=1:num_gratings
            MeanSpike=mean(S(igrat).mean_FRs,2);
            Sort=sort(MeanSpike);
            for i=1:4
            MostActiveNeurons(igrat,i)=find(MeanSpike==Sort(length(Sort)-(i-1)));
            t=0:0.02:1-0.02
            color={[0.4940 0.1840 0.5560],[1 0 0],[0 1 1],[0 0 0]}
            subplot(4,3,igrat)
            plot(t,S(igrat).mean_FRs(MostActiveNeurons(igrat,i),:),'color',color{1,i},'Linewidth',1.5)
            hold all
            b=setdiff([1:size(S(igrat).mean_FRs,1)],MostActiveNeurons(igrat,:));
            a(i)=max(S(igrat).mean_FRs(MostActiveNeurons(igrat,i),:));
            ylim([0 max(a)+0.1])
            end
            for j=1:length(b)
            plot(t,S(igrat).mean_FRs(b(j),:),'color',[0.4660 0.6740 0.1880])
            end
            h=get(gca,'Children');
            Lgnd=legend(h([length(h) length(h)-1 length(h)-2 length(h)-3]),['Neuron ',num2str(NeuronID(MostActiveNeurons(igrat,1)))],['Neuron ',num2str(NeuronID(MostActiveNeurons(igrat,2)))],...
                ['Neuron ',num2str(NeuronID(MostActiveNeurons(igrat,3)))],['Neuron ',num2str(NeuronID(MostActiveNeurons(igrat,4)))],'interpreter', 'latex','Location','northeastoutside');
            xlabel('Time (s)','interpreter', 'latex')
            ylabel('FR(Hz)','interpreter', 'latex');
        end
        clear h
        clear Lgnd
       % fig_name = strcat('FIG_',num2str(imonkey));
       % print(f,fig_name,'-depsc'); 
        clear f 
        MActiveNeurons=reshape(MostActiveNeurons,[1,48]);
        MANeurons=unique(MActiveNeurons)
        for j=1:length(MANeurons)
            counter(j)=length(find(MActiveNeurons==MANeurons(j)))
        end
        counterSort=sort(counter);
        NumNeurons=[];
        for i=1:3
            temp=MANeurons(find(counter==counterSort(length(counterSort)-(i-1))));
            NumNeurons=[NumNeurons temp];
        end
        NumNeurons=unique(NumNeurons);
        NumNeurons=NumNeurons(end-2:end)
        TuningCurve=[];
        Labels=[0:30:330]
        f=figure;
        set(gcf,'units','points','position',[0,0,1000,400])
        for nNeuron=1:length(NumNeurons)
            for igrat=1:num_gratings
                TuningCurve(igrat)=mean(S(igrat).mean_FRs(NumNeurons(nNeuron)));
            end
           subplot(1,3,nNeuron)
           plot(Labels,TuningCurve,'-o','MarkerFaceColor',[0 0 0])
           title(['Neuron ',num2str(NeuronID(NumNeurons(nNeuron)))],'interpreter', 'latex')
           xlim([0 330])
           xlabel('Gratings','interpreter', 'latex');
           ylabel('FR(Hz)','interpreter', 'latex');
        end
        %fig_name = strcat('FIG_',num2str(imonkey+3));
        %print(f,fig_name,'-depsc'); 
        clear S
        clear f
    end
        
   %% Part 2
clear;clc;close all; 

    filenames{1} = 'S_monkey1.mat';
    filenames{2} = 'S_monkey2.mat';
    filenames{3} = 'S_monkey3.mat';
    
    keepNeuronsName{1} = 'keepNeurons1.mat';
    keepNeuronsName{2} = 'keepNeurons2.mat';
    keepNeuronsName{3} = 'keepNeurons3.mat';
    
    dataName{1}='data_monkey1_gratings.mat';
    dataName{2}='data_monkey2_gratings.mat';
    dataName{3}='data_monkey3_gratings.mat';

    monkeys = {'monkey1', 'monkey2', 'monkey3'};
    f=figure;
    set(gcf,'units','points','position',[0,0,600,400])
    for imonkey=1:3
        load(filenames{imonkey});
        load(keepNeuronsName{imonkey});
        load(dataName{imonkey});
        num_gratings=length(S)
        NeuronID=find(keepNeurons==1);
        deletedNeurons=find(keepNeurons==0);
        Channels=data.CHANNELS(:,1);
        MAP=data.MAP;

            for igrat=1:num_gratings
            MeanSpike(:,igrat)=mean(S(igrat).mean_FRs,2);
            end
            gratingMax=zeros(1,length(keepNeurons));
            for j=1:size(MeanSpike,1)
                gratingMax(NeuronID(j))= find(MeanSpike(j,:)==max(MeanSpike(j,:)))
            end     
            MapEdited=reshape(MAP,[100,1]);
            uniqueCH=unique(Channels);
            dif=setdiff(MapEdited,uniqueCH);
            dif(isnan(dif))=[]
            for i=1:length(dif)
            MapEdited(MapEdited==dif(i))=NaN;
            end
            MapEdited=reshape(MapEdited,[10,10]);
            for k=1:length(gratingMax)
                    MapEdited(MapEdited==Channels(k))=gratingMax(k)
            end
             MapEdited(MapEdited==0)=NaN;
             pcolor([MapEdited nan(10,1); nan(1,10+1)]);
             shading flat;
             set(gca, 'ydir', 'reverse');
             c=colorbar
             c.Label.String = 'Gratings';
             title(['Monkey ',num2str(imonkey)],'interpreter', 'latex') 
             
             clear MeanSpike
             clear data
             clear keepNeurons
    end
   % fig_name = strcat('FIG_',num2str(imonkey+12));
   % print(f,fig_name,'-depsc');

      %% Part 3_1
clear;clc;close all; 

    filenames{1} = 'S_monkey1.mat';
    filenames{2} = 'S_monkey2.mat';
    filenames{3} = 'S_monkey3.mat';
    
    keepNeuronsName{1} = 'keepNeurons1.mat';
    keepNeuronsName{2} = 'keepNeurons2.mat';
    keepNeuronsName{3} = 'keepNeurons3.mat';
 
    dataName{1}='data_monkey1_gratings.mat';
    dataName{2}='data_monkey2_gratings.mat';
    dataName{3}='data_monkey3_gratings.mat';

    monkeys = {'monkey1', 'monkey2', 'monkey3'};

    for imonkey=1:3
    load(filenames{imonkey});
    load(keepNeuronsName{imonkey});
    load(dataName{imonkey});
    num_gratings=length(S)
    f=figure;
    for nNeuron=1:length(find(keepNeurons==1))
        TuningCurve=[];
        for igrat=1:num_gratings
            TuningCurve(igrat)=mean(S(igrat).mean_FRs(nNeuron));
        end
            TuningCurveMatrix(:,nNeuron)=TuningCurve;
    end
    r_Sig=corr(TuningCurveMatrix);
   [group1,group2,group3,group4]=rSig_Group(r_Sig);
   
   distance=zeros(length(find(keepNeurons==1)),length(find(keepNeurons==1)));
   for i=1:length(find(keepNeurons==1))-1
       for j=i+1:length(find(keepNeurons==1))
           [xNeuron1,yNeuron1]=find(data.MAP==data.CHANNELS(i));
           [xNeuron2,yNeuron2]=find(data.MAP==data.CHANNELS(j));
           distance(i,j)=(400/1000)*sqrt((yNeuron2-yNeuron1)^2+(xNeuron2-xNeuron1)^2);
        end
   end
   %save(sprintf('distance%d.mat',imonkey),'distance');
   r_sc=zeros(length(find(keepNeurons==1)),length(find(keepNeurons==1)));
   g1=zeros(length(group1),2);
   for k=1:length(group1) 
       Neurons=group1(k,:);
       TuningCurveG=[TuningCurveMatrix(:,Neurons(1)) TuningCurveMatrix(:,Neurons(2))];
       gratMax=[min(find(TuningCurveG(:,1)==max(TuningCurveG(:,1)))) min(find(TuningCurveG(:,2)==max(TuningCurveG(:,2))))];
       count=zeros(2,200);
       for trial=1:200
          count(1,trial)=sum(S(gratMax(1)).trial(trial).counts(Neurons(1),:));
          count(2,trial)=sum(S(gratMax(2)).trial(trial).counts(Neurons(2),:));
       end
       r=corr(count');
       r_sc(Neurons(1),Neurons(2))=r(3);
       g1(k,:)=[distance(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))]
   end
      Mean=0;
   dist1=unique(g1(:,1));
   for m=1:length(dist1)
       index=find(g1(:,1)==dist1(m));
       Sum=0;
       for v=1:length(index)
       Sum=g1(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
   plot(dist1,Mean,'b','Linewidth',1.2);
   
   g2=zeros(length(group2),2);
      for k=1:length(group2) 
       Neurons=group2(k,:);
       TuningCurveG=[TuningCurveMatrix(:,Neurons(1)) TuningCurveMatrix(:,Neurons(2))];
       gratMax=[min(find(TuningCurveG(:,1)==max(TuningCurveG(:,1)))) min(find(TuningCurveG(:,2)==max(TuningCurveG(:,2))))];
       count=zeros(2,200);
       for trial=1:200
          count(1,trial)=sum(S(gratMax(1)).trial(trial).counts(Neurons(1),:));
          count(2,trial)=sum(S(gratMax(2)).trial(trial).counts(Neurons(2),:));
       end
       r=corr(count');
       r_sc(Neurons(1),Neurons(2))=r(3);
       g2(k,:)=[distance(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))];
      end
         Mean=0;
       dist2=unique(g2(:,1));
   for m=1:length(dist2)
       index=find(g2(:,1)==dist2(m));
       Sum=0;
       for v=1:length(index)
       Sum=g2(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
   hold on
   plot(dist2,Mean,'r','Linewidth',1.2);
   
   g3=zeros(length(group3),2);
   for k=1:length(group3) 
       Neurons=group3(k,:);
       TuningCurveG=[TuningCurveMatrix(:,Neurons(1)) TuningCurveMatrix(:,Neurons(2))];
       gratMax=[min(find(TuningCurveG(:,1)==max(TuningCurveG(:,1)))) min(find(TuningCurveG(:,2)==max(TuningCurveG(:,2))))];
       count=zeros(2,200);
       for trial=1:200
          count(1,trial)=sum(S(gratMax(1)).trial(trial).counts(Neurons(1),:));
          count(2,trial)=sum(S(gratMax(2)).trial(trial).counts(Neurons(2),:));
       end
       r=corr(count');
       r_sc(Neurons(1),Neurons(2))=r(3);
       g3(k,:)=[distance(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))];
   end
   Mean=0;
          dist3=unique(g3(:,1));
   for m=1:length(dist3)
       index=find(g3(:,1)==dist3(m));
       Sum=0;
       for v=1:length(index)
       Sum=g3(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
   hold on
   plot(dist3,Mean,'k','Linewidth',1.2);
   
   g4=zeros(length(group4),2);
      for k=1:length(group4) 
       Neurons=group4(k,:);
       TuningCurveG=[TuningCurveMatrix(:,Neurons(1)) TuningCurveMatrix(:,Neurons(2))];
       gratMax=[min(find(TuningCurveG(:,1)==max(TuningCurveG(:,1)))) min(find(TuningCurveG(:,2)==max(TuningCurveG(:,2))))];
       count=zeros(2,200);
       for trial=1:200
          count(1,trial)=sum(S(gratMax(1)).trial(trial).counts(Neurons(1),:));
          count(2,trial)=sum(S(gratMax(2)).trial(trial).counts(Neurons(2),:));
       end
       r=corr(count');
       r_sc(Neurons(1),Neurons(2))=r(3);
       g4(k,:)=[distance(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))];
      end
         Mean=0;
          dist4=unique(g4(:,1));
   for m=1:length(dist4)
       index=find(g4(:,1)==dist4(m));
       Sum=0;
       for v=1:length(index)
       Sum=g4(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
   save(sprintf('rsc_%d.mat',imonkey),'r_sc');
   save(sprintf('r_Sig%d.mat',imonkey),'r_Sig');
   hold on
   plot(dist4,Mean,'g','Linewidth',1.2);
   legend('rsig>0.5','0<rsig<0.5','-0.5<rsig<0','rsig<-0.5','location','northeastoutside');
   title(['Monkey ',num2str(imonkey)],'interpreter', 'latex')
   xlabel('Distance(mm)','interpreter', 'latex');
   ylabel('$r_{sc}$','interpreter', 'latex');
   clear count
   clear gratmax
   clear data
   clear S
    %fig_name = strcat('FIG_',num2str(imonkey+6));
    %print(f,fig_name,'-depsc');
    clear f
    end
    
     %% Part 3_2
    clear;clc;close all;
    rscName{1}='rsc_1.mat';
    rscName{2}='rsc_2.mat';
    rscName{3}='rsc_3.mat';
    
    rsigName{1}='r_Sig1.mat';
    rsigName{2}='r_Sig2.mat';
    rsigName{3}='r_Sig3.mat';
    
    disName{1}='distance1.mat';
    disName{2}='distance2.mat';
    disName{3}='distance3.mat';

  for imonkey=2:2
    load(rscName{imonkey});
    load(rsigName{imonkey});
    load(disName{imonkey});

      group1=0;
      group2=0;
      group3=0;
      group4=0;
      group5=0;
      dis=distance;
      dis(dis>0 & dis<1)=1;
      dis(dis>1 & dis<2)=2;
      dis(dis>2 & dis<3)=3;
      dis(dis>3 & dis<4)=4;
      dis(dis>4 & dis<5)=5;

      for j=1:5
    [N1,N2]=find(dis==j);
    if j==1
    group1=[N1 N2];
    end
    if j==2
        group2=[N1 N2];
    end
     if j==3
        group3=[N1 N2];
     end
     if j==4
        group4=[N1 N2];
     end
     end
        f=figure;
       g1=zeros(length(group1),2);
   for k=1:length(group1) 
       Neurons=group1(k,:);
       g1(k,:)=[r_Sig(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))]
   end
      Mean=0;
   r_sig1=unique(g1(:,1));
   for m=1:length(r_sig1)
       index=find(g1(:,1)==r_sig1(m));
       Sum=0;
       for v=1:length(index)
       Sum=g1(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
      Mean=movmean(Mean,50);

   plot(r_sig1,Mean,'b','Linewidth',1.2);
   
       g2=zeros(length(group2),2);
   for k=1:length(group2) 
       Neurons=group2(k,:);
       g2(k,:)=[r_Sig(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))]
   end
      Mean=0;
   r_sig2=unique(g2(:,1));
   for m=1:length(r_sig2)
       index=find(g2(:,1)==r_sig2(m));
       Sum=0;
       for v=1:length(index)
       Sum=g2(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
      Mean=movmean(Mean,50);
      hold on
   plot(r_sig2,Mean,'r','Linewidth',1.2);

             g3=zeros(length(group3),2);
   for k=1:length(group3) 
       Neurons=group3(k,:);
       g3(k,:)=[r_Sig(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))]
   end
      Mean=0;
   r_sig3=unique(g3(:,1));
   for m=1:length(r_sig3)
       index=find(g3(:,1)==r_sig3(m));
       Sum=0;
       for v=1:length(index)
       Sum=g3(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
      Mean=movmean(Mean,50);

      hold on

   plot(r_sig3,Mean,'k','Linewidth',1.2);

               g4=zeros(length(group4),2);
   for k=1:length(group4) 
       Neurons=group4(k,:);
       g4(k,:)=[r_Sig(Neurons(1),Neurons(2)) r_sc(Neurons(1),Neurons(2))]
   end
      Mean=0;
   r_sig4=unique(g4(:,1));
   for m=1:length(r_sig4)
       index=find(g4(:,1)==r_sig4(m));
       Sum=0;
       for v=1:length(index)
       Sum=g4(index(v),2)+Sum;
       end
       Mean(m)=Sum/length(index);
   end
      Mean=movmean(Mean,50);
      hold on
   plot(r_sig4,Mean,'g','Linewidth',1.2); 
   legend('0<distance<1','1<distance<2','2<distance<3','3<distance<4','location','northeastoutside');
   title(['Monkey ',num2str(imonkey)],'interpreter', 'latex')
   xlabel('$r_{sig}$','interpreter', 'latex');
   ylabel('$r_{sc}$','interpreter', 'latex');
    %fig_name = strcat('FIG_',num2str(imonkey+9));
    %print(f,fig_name,'-depsc');
    clear f
  end

    
  %% Part 3_3
    clear;clc;close all;
    rscName{1}='rsc_1.mat';
    rscName{2}='rsc_2.mat';
    rscName{3}='rsc_3.mat';

    rsigName{1}='r_Sig1.mat';
    rsigName{2}='r_Sig2.mat';
    rsigName{3}='r_Sig3.mat';
    
    disName{1}='distance1.mat';
    disName{2}='distance2.mat';
    disName{3}='distance3.mat';
    for imonkey=1:3
    load(rscName{imonkey});
    load(rsigName{imonkey});
    load(disName{imonkey});
    r_Sig=triu(r_Sig);
    r_Sig(r_Sig==1)=0;
    
      group1=0;group2=0;group3=0;group4=0;group5=0;group6=0;
      group7=0;group8=0;

      dis=distance;
      dis(dis>0 & dis<0.5)=1;
      dis(dis>0.5 & dis<1)=2;
      dis(dis>1 & dis<1.5)=3;
      dis(dis>1.5 & dis<2)=4;
      dis(dis>2 & dis<2.5)=5;
      dis(dis>2.5 & dis<3)=6;
      dis(dis>3 & dis<3.5)=7;
      dis(dis>3.5 & dis<4)=8;

      for j=1:8
    [N1,N2]=find(dis==j);
    if j==1
       group1=[N1 N2];
    end
    if j==2
       group2=[N1 N2];
    end
    if j==3
       group3=[N1 N2];
    end
    if j==4
       group4=[N1 N2];
    end
    if j==5
       group5=[N1 N2];
    end
    if j==6
        group6=[N1 N2];
    end
    if j==7
        group7=[N1 N2];
    end
    if j==8
        group8=[N1 N2];
    end
      end
      Group1=0;Group2=0;Group3=0;Group4=0;Group5=0;Group6=0;
      Group7=0;Group8=0;
      
      r_sig=r_Sig;
      r_sig(r_sig>-1 & r_sig<-0.75)=1;
      r_sig(r_sig>-0.75 & r_sig<-0.5)=2;
      r_sig(r_sig>-0.5 & r_sig<-0.25)=3;
      r_sig(r_sig>-0.25 & r_sig<0)=4;
      r_sig(r_sig>0 & r_sig<0.25)=5;
      r_sig(r_sig>0.25 & r_sig<0.5)=6;
      r_sig(r_sig>0.5 & r_sig<0.75)=7;
      r_sig(r_sig>0.75 & r_sig<1)=8;

    for k=1:8
    [N1,N2]=find(r_sig==k);
    if k==1
       Group1=[N1 N2];
    end
    if k==2
       Group2=[N1 N2];
    end
    if k==3
       Group3=[N1 N2];
    end
    if k==4
       Group4=[N1 N2];
    end
    if k==5
       Group5=[N1 N2];
    end
    if k==6
        Group6=[N1 N2];
    end
    if k==7
        Group7=[N1 N2];
    end
    if k==8
        Group8=[N1 N2];
    end     
    end
    group={group1,group2,group3,group4,group5,group6,group7,group8};
    Group={Group1,Group2,Group3,Group4,Group5,Group6,Group7,Group8};
    for i=1:length(group)
        x=1;
        for j=length(Group):-1:1
            CommGg{x,i}=intersect(group{1,i},Group{1,j},'rows')
            x=x+1;
        end
    end
    for n=1:length(group)
        for m=1:length(Group)
            temp1=CommGg{n,m}
            SumRsc=0;
            for p=1:size(temp1,1)
                temp2=temp1(p,:);
                SumRsc=SumRsc+r_sc(temp2(1),temp2(2));
            end
            MeanRsc(n,m)=SumRsc/size(temp1,1)
        end
    end
    f=figure;
   imagesc([0:0.5:4],[1:-0.25:-1],MeanRsc)
   set(gca,'YDir','normal')
   title(['Monkey ',num2str(imonkey)],'interpreter', 'latex')
   xlabel('Distance(mm)','interpreter', 'latex');
   ylabel('$r_{sig}$','interpreter', 'latex');
   c=colorbar
   c.Label.String = 'r_{sc}';
 %  fig_name = strcat('FIG_',num2str(imonkey+15));
 %  print(f,fig_name,'-depsc');
   clear f;
    end
    
    function [group1,group2,group3,group4]=rSig_Group(r_Sig)
    r_Sig_g=triu(r_Sig);
    r_Sig_g(r_Sig_g==1)=0
    r_Sig_g(r_Sig_g>0.5)=1;
    r_Sig_g(r_Sig_g>0 & r_Sig_g<0.5)=2;
    r_Sig_g(r_Sig_g<0 & r_Sig_g>-0.5)=3;
    r_Sig_g(r_Sig_g<-0.5)=4;
    
    for i=1:4
    [N1,N2]=find(r_Sig_g==i);
    if i==1
    group1=[N1 N2];
    end
    if i==2
        group2=[N1 N2];
    end
     if i==3
        group3=[N1 N2];
     end
     if i==4
        group4=[N1 N2];
    end
    end
    end
      
        
        
        
        
        
        
        
        