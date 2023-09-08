
clear;close all
class = '1h';
fix = load(['E:\',class,'\nor1\',class,'_ct_pointsANDnormal.txt']);
load(['E:\',class,'\nor1\',class,'ct.mat']);
fix(:,7) = normalization_kG';


move = load(['E:\',class,'\nor1\',class,'_mri_pointsANDnormal.txt']);
load(['E:\',class,'\nor1\',class,'mri.mat']);
move(:,7) = normalization_kG';
model = fix';
data = move';
tic;
[RotMat,TransVec,dataOut]=icp7(model,data);
toc;

ct = pointCloud(model(1:3,:)');
mri = pointCloud(dataOut');

figure;pcshowpair(ct,mri,'VerticalAxis','Z','VerticalAxisDir','Down');

pointfix = model';
pointmove = dataOut';

[k1,n]=size(pointfix);
[k2,m]=size(pointmove);
if k1>k2
    pointmove = fix;
    pointfix = dataOut';
    
    [k1,n]=size(pointfix);%k1<k2
    [k2,m]=size(pointmove);
end
datap1=zeros(k2,3);%中间点集
datapp=zeros(k1,3);%对应点集
datapp_index = zeros(k1,1);
distance=zeros(k1,1);
error=zeros(k1,1);%对应点之间误差点集

for i=1:1:k1
    datap1(:,1)=pointmove(:,1)-pointfix(i,1);
    datap1(:,2)=pointmove(:,2)-pointfix(i,2);
    datap1(:,3)=pointmove(:,3)-pointfix(i,3);
    distance=datap1(:,1).^2+datap1(:,2).^2+datap1(:,3).^2;
    [mi,minindex]=min(distance);
    datapp(i,:)=pointmove(minindex,1:3);
    datapp_index(i,1)=minindex;
    error(i)=mi;
end                %用最近距离找对应点
err_mean = mean(error);
err_sum = sum(error);
