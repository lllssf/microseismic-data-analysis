function [R,DH,DM]=getRandSamples(y,features,D1,D2,k)
%---------测试用
%clear;
%[y,features] = signalfeatures();
%k = 8;
%D1 = zeros(0,15);
%D2 = zeros(0,15);
%for t = 1:300
    %if y(t)^2 >= 0.2
      % D1(size(D1,1)+1,:) = features(t,:);
  % else
    %   D2(size(D2,1)+1,:) = features(t,:);
  % end
%end

%---------正式开始
%先产生一个随机数，确定选择的样本R
r = unidrnd(300);
R = features(r,:);
%d1 = zeros(1,size(D1,1));
%d2 = zeros(1,size(D2,1));
for i=1:size(D1,1)
    d1(1,i) = pdist2(R,D1(i,:),'Euclidean');
end
for j=1:size(D2,1)
    d2(1,j) = pdist2(R,D2(j,:),'Euclidean');
end
[v1,L1] = sort(d1);     %d1排序
[v2,L2] = sort(d2);     %d2排序
if y(r)^2 >= 0.2
    H = D1(L1(1,2:k+1),:);      %L1是与R距离最近的编号，赋给H
    M = D2(L2(1,1:k),:);
else
    H = D1(L1(1,1:k),:);
    M = D2(L2(1,2:k+1),:);
end
% 循环计算每2个样本特征之间的特征距离
DH = zeros(k,size(features,2));
DM = zeros(k,size(features,2));
for i = 1:k
    for j = 1:size(features,2)
        DH(i,j) = abs(H(i,j)-R(1,j));
        DM(i,j) = abs(M(i,j)-R(1,j));
    end
end
end