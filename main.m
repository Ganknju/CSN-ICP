%��ά,123ά��xyz���꣬456ά�Ƿ���������

clear;close all
%outside
fix0 = load('E:\code\pointcloud\ct.txt'); 
move0 = load('E:\code\pointcloud\mri.txt');
% fix0 = load('D:\A�о���֮·\project0923\pointcloud\newSide6_ct.txt');
% move0 = load('D:\A�о���֮·\project0923\pointcloud\newSide6_mri.txt');
%normal
fix_normal = load('E:\code\pointcloud\ct_pointsANDnormal.txt'); 
move_normal = load('E:\code\pointcloud\mri_pointsANDnormal.txt');

% tic 
t=clock;
outpoint = fix0;
normal = fix_normal;
fix0_normal = normal_find(outpoint,normal);      %ct�����뷨����ƥ�䣬�����С32387*6
outpoint = move0;
normal = move_normal;
move0_normal = normal_find(outpoint,normal);     %mri�����뷨����ƥ��


%all
fix01 = load('E:\code\pointcloud\ct.txt');
move01 = load('E:\code\pointcloud\mri.txt');

model = fix0_normal';                      %��CT����������Ϊģ�ͣ�ת�ú�����С��Ϊ6*32837
data = move0_normal';                      %��MRI������Ϊ����ƥ�䵽CT�ϣ�ת�ú�����С��Ϊ6*24339
otherdata = move01';                       %otherdata��mri��ȫ���������ݣ�δ��������ƥ�䣩

[RotMat,TransVec,dataOut,otherdataOut]=icp(model,data,otherdata);%icp4()������ת����ƽ�ƾ���
time=etime(clock,t);

%% ����ƥ�����ʾ
%outside��ʾƥ��������ά����ͼ��
ct = pointCloud(model(1:3,:)');
mri = pointCloud(dataOut');
%all��ʾƥ��������ά��������ͼ��
ct_all = pointCloud(fix01);
mri_all = pointCloud(otherdataOut');

% figure;pcshow(mri,'VerticalAxis','Z','VerticalAxisDir','Down');
figure;pcshowpair(mri,ct,'VerticalAxis','Z','VerticalAxisDir','Down');
figure;pcshowpair(mri_all,ct_all,'VerticalAxis','Z','VerticalAxisDir','Down');

%%ȷ���̶������ƶ���ĵ�����ϵ
pointfix = model';
pointmove = dataOut';
[k1,n]=size(pointfix);       % model'=32387*6,k1��model����CT����
[k2,m]=size(pointmove);      %pointmove=32387*3,k2��data��Ӧ����������MRI����
if k1>k2 %��MRI���Ʋ�û��ƥ����CT����
    % pointfix = model';
    pointmove = fix0;
    pointfix = dataOut';%�������ݼ�
    
    [k1,n]=size(pointfix);%k1<k2
    [k2,m]=size(pointmove);
end
datap1=zeros(k2,3);%�м�㼯
datapp=zeros(k1,3);%��Ӧ�㼯
datapp_index = zeros(k1,1);
distance=zeros(k1,1);
error=zeros(k1,1);%��Ӧ��֮�����㼯

for i=1:1:k1
    datap1(:,1)=pointmove(:,1)-pointfix(i,1);
    datap1(:,2)=pointmove(:,2)-pointfix(i,2);
    datap1(:,3)=pointmove(:,3)-pointfix(i,3);    %�ֱ����xyz���ƶ������i���̶���Ĳ�ֵ
    distance=datap1(:,1).^2+datap1(:,2).^2+datap1(:,3).^2;
    [mi,minindex]=min(distance);                 %ȡ������ֵ�������Сֵ
    datapp(i,:)=pointmove(minindex,:);           %��ȡ����������ά��Ϣ
    datapp_index(i,1)=minindex;                  
    error(i)=mi;
end                %����������Ҷ�Ӧ��
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
datap1=zeros(k2,3);%�м�㼯
datapp=zeros(k1,3);%��Ӧ�㼯
datapp_index = zeros(k1,1);
distance=zeros(k1,1);
error=zeros(k1,1);%��Ӧ��֮�����㼯

for i=1:1:k1
    datap1(:,1)=pointmove(:,1)-pointfix(i,1);
    datap1(:,2)=pointmove(:,2)-pointfix(i,2);
    datap1(:,3)=pointmove(:,3)-pointfix(i,3);
    distance=datap1(:,1).^2+datap1(:,2).^2+datap1(:,3).^2;
    [mi,minindex]=min(distance);
    datapp(i,:)=pointmove(minindex,:);
    datapp_index(i,1)=minindex;
    error(i)=mi;
end                %����������Ҷ�Ӧ��
err_mean_all = mean(error);
err_sum_all = sum(error);

