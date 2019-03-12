function [y,newfeatures,feature_index] = featureselect()
%ReliefF特征提取算法
clear;
[y,features] = signalfeatures();
m = 80;     %抽样次数
k = 8;
N = 20;     %运行次数
for i=1:N
    weights(i,:) = ReliefF(y,features,m,k);
end
% 绘制每次计算得到的权重的图像
%for i = 1:N
 %   plot(1:15, weights(i,:));
  %  hold on;
%end
% 计算N次中每个属性的平均值
for i = 1:size(features,2)
    result(1,i) = sum(weights(:,i))/size(weights,1);
end
%plot(result)
%xlabel('number of features');
%ylabel('weight of features');
%grid;

newfeatures = zeros(300,0);
feature_index = zeros(1,0);
for i = 1:size(features,2)
    if result(1,i)>0.05
        newfeatures(:,size(newfeatures,2)+1) = features(:,i);
        feature_index(1,size(feature_index,2)+1) = i;
    end
end
end