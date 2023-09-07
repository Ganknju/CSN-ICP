%六维,123维是xyz坐标，456维是法向量坐标

clear;close all
%outside
fix0 = load('E:\code\pointcloud\ct.txt'); 
move0 = load('E:\code\pointcloud\mri.txt');
% fix0 = load('D:\A研究生之路\project0923\pointcloud\newSide6_ct.txt');
% move0 = load('D:\A研究生之路\project0923\pointcloud\newSide6_mri.txt');
%normal
fix_normal = load('E:\code\pointcloud\ct_pointsANDnormal.txt'); 
move_normal = load('E:\code\pointcloud\mri_pointsANDnormal.txt');

% tic 
t=clock;
outpoint = fix0;
normal = fix_normal;
fix0_normal = normal_find(outpoint,normal);      %ct点云与法向量匹配，矩阵大小32387*6
outpoint = move0;
normal = move_normal;
move0_normal = normal_find(outpoint,normal);     %mri点云与法向量匹配


%all
fix01 = load('E:\code\pointcloud\ct.txt');
move01 = load('E:\code\pointcloud\mri.txt');

model = fix0_normal';                      %以CT点云数据作为模型，转置后矩阵大小变为6*32837
data = move0_normal';                      %以MRI点云作为数据匹配到CT上，转置后矩阵大小变为6*24339
otherdata = move01';                       %otherdata是mri的全部点云数据（未作法向量匹配）

[RotMat,TransVec,dataOut,otherdataOut]=icp(model,data,otherdata);%icp4()返回旋转矩阵、平移矩阵
time=etime(clock,t);

%% 点云匹配后显示
%outside显示匹配后点云三维坐标图像
ct = pointCloud(model(1:3,:)');
mri = pointCloud(dataOut');
%all显示匹配后点云六维坐标数据图像
ct_all = pointCloud(fix01);
mri_all = pointCloud(otherdataOut');

% figure;pcshow(mri,'VerticalAxis','Z','VerticalAxisDir','Down');
figure;pcshowpair(mri,ct,'VerticalAxis','Z','VerticalAxisDir','Down');
figure;pcshowpair(mri_all,ct_all,'VerticalAxis','Z','VerticalAxisDir','Down');

%%确定固定点与移动点的点数关系
pointfix = model';
pointmove = dataOut';
[k1,n]=size(pointfix);       % model'=32387*6,k1是model——CT点数
[k2,m]=size(pointmove);      %pointmove=32387*3,k2是data对应法向量——MRI点数
if k1>k2 %即MRI点云并没有匹配上CT点云
    % pointfix = model';
    pointmove = fix0;
    pointfix = dataOut';%交换数据集
    
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
    datap1(:,3)=pointmove(:,3)-pointfix(i,3);    %分别计算xyz轴移动点与第i个固定点的差值
    distance=datap1(:,1).^2+datap1(:,2).^2+datap1(:,3).^2;
    [mi,minindex]=min(distance);                 %取遍历差值距离的最小值
    datapp(i,:)=pointmove(minindex,:);           %提取出这个点的六维信息
    datapp_index(i,1)=minindex;                  
    error(i)=mi;
end                %用最近距离找对应点
err_mean = mean(error);
err_sum = sum(error);

%all
pointfix = fix01;
pointmove = otherdataOut';

[k1,n]=size(pointfix);
[k2,m]=size(pointmove);
if k1>k2
   pointmove = fix01;
pointfix = otherdataOut';

[k1,n]=size(pointfix);
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
    datapp(i,:)=pointmove(minindex,:);
    datapp_index(i,1)=minindex;
    error(i)=mi;
end                %用最近距离找对应点
err_mean_all = mean(error);
err_sum_all = sum(error);

