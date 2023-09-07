function [TR,TT,data,otherdata] = icp(model,data,otherdata,maxIter,minIter,critFun,thres)

data_nor = data(4:6,:);
model_nor = model(4:6,:);
data_xor = data(1:3,:);
model_xor = model(1:3,:);

data(4:6,:) = [];
model(4:6,:) = [];            %删除data与model矩阵的4到6行
% Check input arguments

if nargin<3                     %nargin可返回函数输入参数的数量
   
    error('To few input arguments');
    
elseif nargin<7
    
    thres=2e-5;                     % threshold to stop icp iterations
    if nargin<6
        critFun=0;                  % critFun method, LS
        if nargin<5
            minIter=5;              % min number of icp iterations
            if nargin<4
                maxIter=100;        % max number of icp iterations
            end
        end
    end
    
end

if or(isempty(model),isempty(data))
    error('Something is wrong with the model points and data points');
end

% Use default values            %使用默认值

if isempty(maxIter)             %本代码使用时仅输入前三个参数，即后三个参数均为空字符串，则需要使用默认值
    maxIter=100;
end

if isempty(minIter)
    minIter=5;
end

if isempty(critFun)
    critFun=0;
end

if isempty(thres)
    thres=2e-5;
end

% Size of model points and data points                 %比较数据集大小

if (size(model,2)<size(model,1))                              %size(X,1),返回矩阵X的行数；size(X,2),返回矩阵X的列数；model=6*32837
    mTranspose=true;
    m=size(model,2);
    M=size(model,1);
else
    mTranspose=false;
    m=size(model,1);                                                %正常情况下，m＝6，M=32837
    M=size(model,2);
end

if (size(data,2)<size(data,1))                                  %data=6x24339
    data=data';
end

if m~=size(data,1)
    error('The dimension of the model points and data points must be equal');
end

N=size(data,2);%m是model的行数即维数，M是model的列数即数据点数，N是data的列数即数据点数

% Create closest point search structure       %最近点迭代

 if m<4
    if mTranspose
        DT=delaunayTriangulation(model);
    else
        DT=delaunayTriangulation(model');       %DT三角网格划分
    end
else
    DT=[];
    resid=zeros(N,1);
    vi=ones(N,1);
end

% Initiate weights (Only for robust criterion) 初始化权重

if critFun>0
    wghs=ones(N,1);
end

% Initiate transformation

TR=eye(m);      %旋转矩阵设置为m维单位矩阵
TT=zeros(m,1);  %平移矩阵设置为m*1的列向量

% Start the ICP algorithm

res=9e99;                                                                 %9e99＝9*10^99
res_angel = 9e99;

tree = KDTreeSearcher(model');  % Kd树算法通过将K维空间中的n个点递归地分割为二叉树来对n×K数据集进行分区。
radius = 0.05;
for iter=1:maxIter                                                     %迭代最大值为100
    
     oldres=res;
    
    % Find closest model points to data points
    
    [vi,~] = nearestNeighbor(DT,data');  %DT三角划分 ；data'查询点    data移动匹配模型model
    %vi:最接近查询点的DT顶点ID
    %resid：对应的欧几里得距离
%     %%法线去除极大偏差
%     nor_bias = data_nor-model_nor(:,vi);    
%     nor_bias_res =  nor_bias(1,:).^2 + nor_bias(2,:).^2 + nor_bias(3,:).^2;
    
    % Find neighbors within the fixed radius
    idx = rangesearch(tree,model_xor(:,vi)', radius);      %寻找tree集合中与query集合中的点的距离小于radius的点，即最近点
    idxs = idx{1};                               %idx{1}代表元胞数组idx第一行的元胞
    %%利用坐标点和邻域点
    neighbors = [idxs(:,1) idxs(:,2) idxs(:,3)];
    %%法向量内积
    nor_bias_res = model_nor(:,vi)'.*neighbors;

    vi1 = vi;
    N1=N;
%     vi1(nor_bias_index)=[];
%     data1(:,nor_bias_index)=[];
%     resid(nor_bias_index)=[];
%     N1=size(data1,2);
    % Find transformation
            res=mean(nor_bias_res);
            
            med=mean(data,2);   %[3,1] data的xyz三个方向上均值,mean(a,2)对每行求均值
            if mTranspose  %pass
                mem=mean(model(vi1,:),1);
                C=data*model(vi1,:)-(N1*med)*mem;
                [U,~,V]=svd(C);
                Ri=V*U';
                
                d=eig(C);  %求矩阵C的全部特征值,构成列向量d
                D=sort(d,'descend'); %以降序排序 
                
                
                if det(Ri)<0
                    V(:,end)=-V(:,end);
                    Ri=V*U';
                end
                Ti=mem'-Ri*med;
            else
                mem=mean(model(:,vi1),2);  %根据data查询点得到的model点集 在xyz三个方向上的均值
                C=data*model(:,vi1)'-(N1*med)*mem';  %data*model（data对应的）-N*两者的xyz均值，N是data的数据 
                %%对称矩阵aaT进行svd奇异值分解
                [U,~,V]=svd(C);
                Ri=V*U';
                if det(Ri)<0               %det（）求行列式
                    V(:,end)=-V(:,end);    
                    Ri=V*U';
                end
                Ti=mem-Ri*med;             %model点集在xyz三个方向上的均值-Ri*data的xyz三个方向上均值
            end

    data=Ri*data;                       % Apply transformation 就是重复迭代！
    for i=1:m
        data(i,:)=data(i,:)+Ti(i);      %m是model列数
    end
    otherdata=Ri*otherdata;                       % Apply transformation
    for i=1:m
        otherdata(i,:)=otherdata(i,:)+Ti(i);      %
    end
    
    TR=Ri*TR;                           % Update transformation
    TT=Ri*TT+Ti;                        %
    
    if iter >= minIter
        if abs(oldres-res) < thres
            break
        end
    end
end

