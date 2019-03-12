function weights = ReliefF(y,features,m,k)
% Input： k--最邻近样本个数
%        m--抽样次数
% Output： weights--特征权重

%-----测试
%clear;
%[y,features] = signalfeatures();
%m=80;
%k=8;

%------------先将数据集分为两类，可以加快计算速度
D1 = zeros(0,size(features,2));
D2 = zeros(0,size(features,2));
for t = 1:300
    if y(t)^2 >= 0.2
        D1(size(D1,1)+1,:) = features(t,:);
    else
        D2(size(D2,1)+1,:) = features(t,:);
    end
end
weights = zeros(1,size(features,2));
%进行m次循环选择操作
for i = 1:m 
    %随机选择一个样本R
    [R,DH,DM] = getRandSamples(y,features,D1,D2,k);
    %更新特征权重值
    for j = 1:length(weights)
        weights(1,j) = weights(1,j)-sum(DH(:,j))/(k*m)+sum(DM(:,j))/(k*m);
    end
end
weights = abs(weights);
end