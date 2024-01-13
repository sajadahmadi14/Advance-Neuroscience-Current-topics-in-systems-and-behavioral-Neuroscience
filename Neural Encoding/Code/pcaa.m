
clear all; close all; clc;
load UnitsData.mat;

for n = 1:481
    [h(n,:),h1(n,:),h2(n,:),h3(n,:),h4(n,:),...
     h5(n,:),h6(n,:)]=psth(n,Unit);
end

h=reshape(h,481,53)';
h=h(2:52,:);
[~,pcc,var11]=pca(h);

h1=reshape(h1,481,53)';
h1=h1(2:52,:);
[~,pcc1,var1]=pca(h1);
h2=reshape(h2,481,53)';
h2=h2(2:52,:);
[~,pcc2,var2]=pca(h2);
h3=reshape(h3,481,53)';
h3=h3(2:52,:);
[~,pcc3,var3]=pca(h3);
h4=reshape(h4,481,53)';
h4=h4(2:52,:);
[~,pcc4,var4]=pca(h4);
h5=reshape(h5,481,53)';
h5=h5(2:52,:);
[~,pcc5,var5]=pca(h5);
h6=reshape(h6,481,53)';
h6=h6(2:52,:);
[~,pcc6,var6]=pca(h6);
%hold on
var=[var1 var2 var3 var4 var5 var6];
figure;
plot(pcc(:,1),pcc(:,2));
xlabel('PC1'); ylabel('PC2');title('ALL');
figure;
subplot(2,3,1);
plot(pcc1(:,1),pcc1(:,2));
xlabel('PC1'); ylabel('PC2');title('Condition1');
subplot(2,3,2);
plot(pcc2(:,1),pcc2(:,2));
xlabel('PC1'); ylabel('PC2');title('Condition2');
subplot(2,3,3);
plot(pcc3(:,1),pcc3(:,2));
xlabel('PC1'); ylabel('PC2');title('Condition3');
subplot(2,3,4);
plot(pcc4(:,1),pcc4(:,2));
xlabel('PC1'); ylabel('PC2');title('Condition4');
subplot(2,3,5);
plot(pcc5(:,1),pcc5(:,2));
xlabel('PC1'); ylabel('PC2');title('Condition5');
subplot(2,3,6);
plot(pcc6(:,1),pcc6(:,2));
xlabel('PC1'); ylabel('PC2');title('Condition6');



figure;
plot3(pcc(:,1),pcc(:,2),pcc(:,3));
xlabel('PC1'); ylabel('PC2');zlabel('PC3');title('ALL');


figure;
subplot(2,3,1);
plot3(pcc1(:,1),pcc1(:,2),pcc1(:,3));
xlabel('PC1'); ylabel('PC2');zlabel('PC3');title('Condition1');
subplot(2,3,2);
plot3(pcc2(:,1),pcc2(:,2),pcc2(:,3));
xlabel('PC1'); ylabel('PC2');zlabel('PC3');title('Condition2');
subplot(2,3,3);
plot3(pcc3(:,1),pcc3(:,2),pcc3(:,3));
xlabel('PC1'); ylabel('PC2');zlabel('PC3');title('Condition3');
subplot(2,3,4);
plot3(pcc4(:,1),pcc4(:,2),pcc4(:,3));
xlabel('PC1'); ylabel('PC2');zlabel('PC3');title('Condition4');
subplot(2,3,5);
plot3(pcc5(:,1),pcc5(:,2),pcc5(:,3));
xlabel('PC1'); ylabel('PC2');zlabel('PC3');title('Condition5');
subplot(2,3,6);
plot3(pcc6(:,1),pcc6(:,2),pcc6(:,3));
xlabel('PC1'); ylabel('PC2');zlabel('PC3');title('Condition6');


figure;
hist(var11,50);
title('variance of PCs','ALL');xlabel('PCs');


figure;
subplot(2,3,1);
hist(var1,50);
title('variance of PCs','Condition1');xlabel('PCs');
subplot(2,3,2);
hist(var2,50);
title('variance of PCs','Condition2');xlabel('PCs');
subplot(2,3,3);
hist(var3,50);
title('variance of PCs','Condition3');xlabel('PCs');
subplot(2,3,4);
hist(var4,50);
title('variance of PCs','Condition4');xlabel('PCs');
subplot(2,3,5);
hist(var5,50);
title('variance of PCs','Condition5');xlabel('PCs');
subplot(2,3,6);
hist(var6,50);
title('variance of PCs','Condition6');xlabel('PCs');