function [TR,TT,data,otherdata] = icp(model,data,otherdata,maxIter,minIter,critFun,thres)

data_nor = data(4:6,:);
model_nor = model(4:6,:);
data_xor = data(1:3,:);
model_xor = model(1:3,:);

data(4:6,:) = [];
model(4:6,:) = [];            %ɾ��data��model�����4��6��
% Check input arguments

if nargin<3                     %nargin�ɷ��غ����������������
   
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

% Use default values            %ʹ��Ĭ��ֵ

if isempty(maxIter)             %������ʹ��ʱ������ǰ������������������������Ϊ���ַ���������Ҫʹ��Ĭ��ֵ
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

% Size of model points and data points                 %�Ƚ����ݼ���С

if (size(model,2)<size(model,1))                              %size(X,1),���ؾ���X��������size(X,2),���ؾ���X��������model=6*32837
    mTranspose=true;
    m=size(model,2);
    M=size(model,1);
else
    mTranspose=false;
    m=size(model,1);                                                %��������£�m��6��M=32837
    M=size(model,2);
end

if (size(data,2)<size(data,1))                                  %data=6x24339
    data=data';
end

if m~=size(data,1)
    error('The dimension of the model points and data points must be equal');
end

N=size(data,2);%m��model��������ά����M��model�����������ݵ�����N��data�����������ݵ���

% Create closest point search structure       %��������

 if m<4
    if mTranspose
        DT=delaunayTriangulation(model);
    else
        DT=delaunayTriangulation(model');       %DT�������񻮷�
    end
else
    DT=[];
    resid=zeros(N,1);
    vi=ones(N,1);
end

% Initiate weights (Only for robust criterion) ��ʼ��Ȩ��

if critFun>0
    wghs=ones(N,1);
end

% Initiate transformation

TR=eye(m);      %��ת��������Ϊmά��λ����
TT=zeros(m,1);  %ƽ�ƾ�������Ϊm*1��������

% Start the ICP algorithm

res=9e99;                                                                 %9e99��9*10^99
res_angel = 9e99;

tree = KDTreeSearcher(model');  % Kd���㷨ͨ����Kά�ռ��е�n����ݹ�طָ�Ϊ����������n��K���ݼ����з�����
radius = 0.05;
for iter=1:maxIter                                                     %�������ֵΪ100
    
     oldres=res;
    
    % Find closest model points to data points
    
    [vi,~] = nearestNeighbor(DT,data');  %DT���ǻ��� ��data'��ѯ��    data�ƶ�ƥ��ģ��model
    %vi:��ӽ���ѯ���DT����ID
    %resid����Ӧ��ŷ����þ���
%     %%����ȥ������ƫ��
%     nor_bias = data_nor-model_nor(:,vi);    
%     nor_bias_res =  nor_bias(1,:).^2 + nor_bias(2,:).^2 + nor_bias(3,:).^2;
    
    % Find neighbors within the fixed radius
    idx = rangesearch(tree,model_xor(:,vi)', radius);      %Ѱ��tree��������query�����еĵ�ľ���С��radius�ĵ㣬�������
    idxs = idx{1};                               %idx{1}����Ԫ������idx��һ�е�Ԫ��
    %%���������������
    neighbors = [idxs(:,1) idxs(:,2) idxs(:,3)];
    %%�������ڻ�
    nor_bias_res = model_nor(:,vi)'.*neighbors;

    vi1 = vi;
    N1=N;
%     vi1(nor_bias_index)=[];
%     data1(:,nor_bias_index)=[];
%     resid(nor_bias_index)=[];
%     N1=size(data1,2);
    % Find transformation
            res=mean(nor_bias_res);
            
            med=mean(data,2);   %[3,1] data��xyz���������Ͼ�ֵ,mean(a,2)��ÿ�����ֵ
            if mTranspose  %pass
                mem=mean(model(vi1,:),1);
                C=data*model(vi1,:)-(N1*med)*mem;
                [U,~,V]=svd(C);
                Ri=V*U';
                
                d=eig(C);  %�����C��ȫ������ֵ,����������d
                D=sort(d,'descend'); %�Խ������� 
                
                
                if det(Ri)<0
                    V(:,end)=-V(:,end);
                    Ri=V*U';
                end
                Ti=mem'-Ri*med;
            else
                mem=mean(model(:,vi1),2);  %����data��ѯ��õ���model�㼯 ��xyz���������ϵľ�ֵ
                C=data*model(:,vi1)'-(N1*med)*mem';  %data*model��data��Ӧ�ģ�-N*���ߵ�xyz��ֵ��N��data������ 
                %%�Գƾ���aaT����svd����ֵ�ֽ�
                [U,~,V]=svd(C);
                Ri=V*U';
                if det(Ri)<0               %det����������ʽ
                    V(:,end)=-V(:,end);    
                    Ri=V*U';
                end
                Ti=mem-Ri*med;             %model�㼯��xyz���������ϵľ�ֵ-Ri*data��xyz���������Ͼ�ֵ
            end

    data=Ri*data;                       % Apply transformation �����ظ�������
    for i=1:m
        data(i,:)=data(i,:)+Ti(i);      %m��model����
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

