function [TR,TT,data] = icp7(model,data,maxIter,minIter,critFun,thres)

%%
euc_ther = 3/4;
nor_ther = 1/2;

%%
data_nor = data(4:6,:);
model_nor = model(4:6,:);
data_cur = data(7,:);
model_cur = model(7,:);

[azi_data,ele_data,r_data] = cart2sph(data_nor(1,:),data_nor(2,:),data_nor(3,:));  %azimuth,elevation,r分别表示方位角、仰角、和半径
[azi_model,ele_model,r_model] = cart2sph(model_nor(1,:),model_nor(2,:),model_nor(3,:));
r_data = data_cur;
r_model = model_cur;
data(4:7,:) = [];
model(4:7,:) = [];
% Check input arguments

if nargin<2
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

% Use default values
if isempty(maxIter)
    maxIter=200;
end

if isempty(minIter)
    minIter=5;
end

if isempty(critFun)
    critFun=0;
end

if isempty(thres)
    thres=1e-5;
end

% Size of model points and data points

if (size(model,2)<size(model,1))
    mTranspose=true;
    m=size(model,2);
    M=size(model,1);
else
    mTranspose=false;
    m=size(model,1);
    M=size(model,2);
end

if (size(data,2)<size(data,1))
    data=data';
end

if m~=size(data,1)
    error('The dimension of the model points and data points must be equal');
end

N=size(data,2);

% Create closest point search structure

 if m<4
    if mTranspose
        DT=delaunayTriangulation(model);
    else
        DT=delaunayTriangulation(model'); %DT三角网格划分
    end
else
    DT=[];
    resid=zeros(N,1);
    vi=ones(N,1);
end

% Initiate weights (Only for robust criterion)

if critFun>0
    wghs=ones(N,1);
end

% Initiate transformation

TR=eye(m);  %单位矩阵
TT=zeros(m,1);

% Start the ICP algorithm

res=9e99;
res_angel = 9e99;

for iter=1:maxIter
    
    oldres=res;
    
    % Find closest model points to data points
    
    [vi0,resid] = nearestNeighbor(DT,data');  %DT三角划分 ；data'查询点
    %vi:最接近查询点的DT顶点ID
    %resid：对应的欧几里得距离
    %法线去除极大偏差
    radius = 0.015;%0.015
    normal_base = data_nor';
    tree_source = model';
    tree_normal = model_nor';
    [vi] = find_neigh_normal(normal_base,vi0,tree_source,tree_normal,radius);
    
    A1 = azi_data;
    A2 = ele_data;
    A3 = r_data;
    B1 = azi_model(:,vi);
    B2 = ele_model(:,vi);
    B3 = r_model(:,vi);
    nor_bias_res = (A1-B1)+(A2-B2);  
    min_norbias = min(nor_bias_res);
    max_norbias = max(nor_bias_res);
    nor_bias_index = find(nor_bias_res>(min_norbias+2/3*(max_norbias-min_norbias)));
    
    %if A3~=Null && B3~=Null
        r_bias_res = A3-B3;
        min_rbias = min(r_bias_res);
        max_rbias = max(r_bias_res);
        r_bias_index = find(r_bias_res>(min_rbias+2/3*(max_rbias-min_rbias)));
    %end
    
    %欧式距离去除极大偏差
    euc_bias = data-model(:,vi);
    euc_bias_res =  euc_bias(1,:).^2 + euc_bias(2,:).^2 + euc_bias(3,:).^2;
    min_eucbias = min(euc_bias_res);
    max_eucbias = max(euc_bias_res);
    euc_bias_index = find(euc_bias_res>(min_eucbias+2/3*(max_eucbias-min_eucbias)));

    vi1 = vi;
    data1 = data;

    bias_index = zeros(1,length(euc_bias_index)+length(nor_bias_index)+length(r_bias_index));
    bias_index(1,1:length(euc_bias_index)) = euc_bias_index;
    bias_index(1,length(euc_bias_index)+1:length(euc_bias_index)+length(nor_bias_index)) = nor_bias_index;
    bias_index(1,length(euc_bias_index)+length(nor_bias_index)+1:length(euc_bias_index)+length(nor_bias_index)+length(r_bias_index)) = r_bias_index;
    bias_index = unique(bias_index);
    vi1(bias_index)=[];
    data1(:,bias_index)=[];
    resid(bias_index)=[];
    N1=size(data1,2);



    %% Find transformation
    res=mean(resid.^2)+mean(nor_bias_res)+mean(r_bias_res);  %res：距离的平方和取均值（一个值）
    med=mean(data1,2);   %[3,1] data的xyz三个方向上均值
    if mTranspose  %pass
        mem=mean(model(vi1,:),1);
        C=data1*model(vi1,:)-(N1*med)*mem;
        [U,~,V]=svd(C);
        Ri=V*U';
        if det(Ri)<0
            V(:,end)=-V(:,end);
            Ri=V*U';
        end
        Ti=mem'-Ri*med;
    else
        mem=mean(model(:,vi1),2);  %根据data查询点得到的model点集 在xyz三个方向上的均值
        C=data1*model(:,vi1)'-(N1*med)*mem';  %data*model（data对应的）-N*两者的xyz均值
        [U,~,V]=svd(C);
        Ri=V*U';
        if det(Ri)<0
            V(:,end)=-V(:,end);
            Ri=V*U';
        end
        Ti=mem-Ri*med;
    end

    data=Ri*data;                       % Apply transformation
    for i=1:m
        data(i,:)=data(i,:)+Ti(i);      %
    end
    
    TR=Ri*TR;                           % Update transformation
    TT=Ri*TT+Ti;                        %
    
    if iter >= minIter
        if abs(oldres-res) < thres
            break
        end
    end
    
end

