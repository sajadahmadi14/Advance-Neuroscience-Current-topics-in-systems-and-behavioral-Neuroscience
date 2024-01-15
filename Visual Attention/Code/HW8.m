%% Sajad AhmadiNabi _ ID: 400206584
%   Assignment_8



%% Part 1-a 
close all;clc;clear;

datafolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\DATA\hp';
stimfolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI';

ImageNumbers=randi(1003,[1 6]);
figure;
set(gcf,'units','points','position',[0,0,1600,650])
for i=1:length(ImageNumbers)
    subplot(2,3,i)
    showEyeData(datafolder, stimfolder ,ImageNumbers(i))
    title(sprintf("Image Number : %d",ImageNumbers(i)),'interpreter','latex')
end
print -depsc part1.eps

%% Part 1_b
close all;clc;clear;

datafolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\DATA\hp';
stimfolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI';

showEyeDataAcrossUsers(stimfolder , 6)

%% Part 2 
close all;clc;clear;

datafolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\DATA\hp';
stimfolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI';

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\FaceDetect');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel\horizon code');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel\FelzenszwalbDetectors');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\LabelMeToolbox-master\LabelMeToolbox-master\features');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\LabelMeToolbox-master\LabelMeToolbox-master\imagemanipulation');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\matlabPyrTools-master\matlabPyrTools-master');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\matlabPyrTools-master\matlabPyrTools-master\MEX');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\SaliencyToolbox2.3\SaliencyToolbox');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\voc-release-3.1-win-master\voc-release-3.1-win-master');

files=dir(fullfile(stimfolder,'*.jpeg'));
[filenames{1:size(files,1)}] = deal(files.name);
randomNum=randi(1003,[1 6]);
featureChar={'Subband','Itti','Color','Torralba','Horizon','Object','DistToCenter','All'};
for j=1:length(randomNum)
imgFile =filenames{1,randomNum(j)}  ;
close all;
f=figure;
set(gcf,'units','points','position',[0,0,1600,650])
subplot(3,3,1)
img = imread(imgFile);
imshow(imgFile);
title('Original Image','interpreter','latex');
for i=1:length(featureChar)
    subplot(3,3,i+1)
    saliencyMap =saliency(imgFile,i);
    title(sprintf("Feature: %s ",featureChar{1,i}),'interpreter','latex')
end
    fig_name = strcat('Part23_',num2str(j));
    print(f,fig_name,'-depsc'); 
end

%% Part 3_a
close all;clc;clear;

stimfolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI';

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\FaceDetect');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel\horizon code');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel\FelzenszwalbDetectors');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\LabelMeToolbox-master\LabelMeToolbox-master\features');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\LabelMeToolbox-master\LabelMeToolbox-master\imagemanipulation');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\matlabPyrTools-master\matlabPyrTools-master');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\matlabPyrTools-master\matlabPyrTools-master\MEX');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\SaliencyToolbox2.3\SaliencyToolbox');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\voc-release-3.1-win-master\voc-release-3.1-win-master');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\gbvs\gbvs\util');

files=dir(fullfile(stimfolder,'*.jpeg'));
[filenames{1:size(files,1)}] = deal(files.name);
randomNum=randi(1003,[1 4]);
featureChar={'Subband','Itti','Color','Torralba','Horizon','Object','DistToCenter','All'};
subjects={'ajs','CNG','emb','ems','ff','hp','jcw','jw','kae','krl','po','tmj','tu','ya','zb'};
for k=1:length(randomNum)
imgFile =filenames{1,randomNum(k)};
imgStruct{k} = imread(imgFile);
for i=1:length(featureChar)
    saliencyMap{i,k} =saliency(imgFile,i);
    close 
end
% load fixation maps of any subject
filename=filenames{1,randomNum(k)}(1:end-5);
for j=1:length(subjects)
  datafolder = strcat('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\DATA\',subjects{1,j});
  SubData{k,j} = load(fullfile(datafolder, filename));
  fns=fieldnames(SubData{k,j});
  [SaccData{k,j},~,~]=checkFixations(SubData{k,j}.(fns{1}).DATA.eyeData);
end
for j=1:length(subjects)
for i=1:length(featureChar)
    ROC_first.sub{k,j}.feature{i}=rocScoreSaliencyVsFixations(saliencyMap{i,k}, rmmissing(SaccData{k,j}(1:floor(length(SaccData{k,j})/2),1)), rmmissing(SaccData{k,j}(1:floor(length(SaccData{k,j})/2),2)), size(imgStruct{k}));
    ROC_second.sub{k,j}.feature{i}=rocScoreSaliencyVsFixations(saliencyMap{i,k}, rmmissing(SaccData{k,j}(floor(length(SaccData{k,j})/2)+1:length(SaccData{k,j}),1)), rmmissing(SaccData{k,j}(floor(length(SaccData{k,j})/2)+1:length(SaccData{k,j}),2)), size(imgStruct{k}));
end
end
end

figure;
set(gcf,'units','points','position',[0,0,1600,650])
for k=1:3
    if k==1
        v=1;
    else if k==2
        v=4;
        else if k==3
                v=7;
            end
        end
    end
subplot(3,3,v)
imshow(imgStruct{k})
title('Original Image','interpreter','latex')
subplot(3,3,v+1)
imshow(saliencyMap{8,k})
hold on 
scatter(rmmissing(SaccData{k,6}(1:floor(length(SaccData{k,6})/2),1)),rmmissing(SaccData{k,6}(1:floor(length(SaccData{k,6})/2),2)),'filled','r')
title('First 1.5s','interpreter','latex')
subplot(3,3,v+2)
imshow(saliencyMap{8,k})
hold on 
scatter(rmmissing(SaccData{k,6}(floor(length(SaccData{k,6})/2)+1:length(SaccData{k,6}),1)),rmmissing(SaccData{k,6}(floor(length(SaccData{k,6})/2)+1:length(SaccData{k,6}),2)),'filled','r')
title('Second 1.5s','interpreter','latex')
end
print -depsc part3_a.eps

for imageNum=1:3
MeanRocFirst=zeros(1,8);
MeanRocSecond=zeros(1,8);
for sub=1:15
MeanRocFirst=MeanRocFirst+cell2mat(ROC_first.sub{imageNum,sub}.feature);
MeanRocSecond=MeanRocSecond+cell2mat(ROC_second.sub{imageNum,sub}.feature);
end
MeanRocFirst_img{imageNum}=MeanRocFirst/15;
MeanRocSecond_img{imageNum}=MeanRocSecond/15;
end
%% Part 3_b
close all;clc;clear;

stimfolder = 'J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI';

addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\FaceDetect');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel\horizon code');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\JuddSaliencyModel\JuddSaliencyModel\FelzenszwalbDetectors');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\LabelMeToolbox-master\LabelMeToolbox-master\features');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\LabelMeToolbox-master\LabelMeToolbox-master\imagemanipulation');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\matlabPyrTools-master\matlabPyrTools-master');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\matlabPyrTools-master\matlabPyrTools-master\MEX');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\SaliencyToolbox2.3\SaliencyToolbox');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\ALLSTIMULI');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\voc-release-3.1-win-master\voc-release-3.1-win-master');
addpath('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\gbvs\gbvs\util');

files=dir(fullfile(stimfolder,'*.jpeg'));
[filenames{1:size(files,1)}] = deal(files.name);
randomNum=randi(1003,[1 100]);
featureChar={'Subband','Itti','Color','Torralba','Horizon','Object','DistToCenter','All'};
subjects={'ajs','CNG','emb','ems','ff','hp','jcw','jw','kae','krl','po','tmj','tu','ya','zb'};
for k=1:length(randomNum)
imgFile =filenames{1,randomNum(k)};
imgStruct{k} = imread(imgFile);
for i=1:length(featureChar)
    saliencyMap{i,k} =saliency(imgFile,i);
    close 
end
filename=filenames{1,randomNum(k)}(1:end-5);
for j=1:length(subjects)
  datafolder = strcat('J:\neuroengineering\semester 2\Advance Neuroscience\HW\HW08\Data\Eye tracking database\Eye tracking database\DATA\',subjects{1,j});
  SubData{k,j} = load(fullfile(datafolder, filename));
  fns=fieldnames(SubData{k,j});
  [SaccData{k,j},~,~]=checkFixations(SubData{k,j}.(fns{1}).DATA.eyeData);
end
for j=1:length(subjects)
for i=1:length(featureChar)
    ROC_first.sub{k,j}.feature{i}=rocScoreSaliencyVsFixations(saliencyMap{i,k}, rmmissing(SaccData{k,j}(1:floor(length(SaccData{k,j})/2),1)), rmmissing(SaccData{k,j}(1:floor(length(SaccData{k,j})/2),2)), size(imgStruct{k}));
    ROC_second.sub{k,j}.feature{i}=rocScoreSaliencyVsFixations(saliencyMap{i,k}, rmmissing(SaccData{k,j}(floor(length(SaccData{k,j})/2)+1:length(SaccData{k,j}),1)), rmmissing(SaccData{k,j}(floor(length(SaccData{k,j})/2)+1:length(SaccData{k,j}),2)), size(imgStruct{k}));
end
end
end

for imageNum=1:length(randomNum)
MeanRocFirst=zeros(1,8);
MeanRocSecond=zeros(1,8);
for sub=1:15
MeanRocFirst=MeanRocFirst+cell2mat(ROC_first.sub{imageNum,sub}.feature);
MeanRocSecond=MeanRocSecond+cell2mat(ROC_second.sub{imageNum,sub}.feature);
end
MeanRocFirst_img{imageNum}=MeanRocFirst/15;
MeanRocSecond_img{imageNum}=MeanRocSecond/15;
end

for imageNum=1:length(randomNum)
    diffRoc_img{1,imageNum}=MeanRocSecond_img{1,imageNum}-MeanRocFirst_img{1,imageNum};
    for i=1:8
        if diffRoc_img{1,imageNum}(1,i)<0
           diffRoc_img{1,imageNum}(1,i)=0;
        else if diffRoc_img{1,imageNum}(1,i)>0
                diffRoc_img{1,imageNum}(1,i)=1;
            end
        end
    end
end

perRoc=zeros(1,8);
for imageNum=1:length(randomNum)
   perRoc=perRoc+diffRoc_img{1,imageNum};
end

RocMeanFirst=zeros(1,8);
RocMeanSecond=zeros(1,8);
for imageNum=1:length(randomNum)
   RocMeanFrist=RocMeanFirst+MeanRocFirst_img{1,imageNum};
   RocMeanSecond=RocMeanSecond+MeanRocSecond_img{1,imageNum};
end
   RocMeanFist=RocMeanFist/100;
   RocMeanSecond=RocMeanSecond/100;

    




