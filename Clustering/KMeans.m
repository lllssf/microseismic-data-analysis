function remark = KMeans(features,N)
% Input: N -- numbers of classes
% Output: remark -- 带分类标号的数据
[m, n] = size(features);     %m是数据个数，n是数据维数
ma = zeros(n);      %每一维最大的数
mi = zeros(n);      %每一维最小的数
u = zeros(N,n);     %随机初始化，最终迭代到每一类的中心位置
for i = 1:n
    ma(i) = max(features(:,i));
    mi(i) = min(features(:,i));
    for j = 1:N
        u(j,i) = ma(i)+(mi(i)-ma(i))*rand();    %随机初始化
    end
end

while 1
    pre_u = u;      %上一次求得的中心位置
    for i = 1:N
        tmp{i} = [];
        for j = 1:m
            tmp{i} = [tmp{i};features(j,:)-u(i,:)];
        end
    end
    quan = zeros(m,N);
    for i=1:m
        c=[];
        for j = 1:N
            c = [c norm(tmp{j}(i,:))];
        end
        [junk, index]=min(c);
        quan(i,index) = norm(tmp{index}(i,:));
    end
    
    for i = 1:N
        for j=1:n
            u(i,j) = sum(quan(:,i).*features(:,j))/sum(quan(:,i));
        end
    end
    
    if norm(pre_u-u)<0.1
        break;
    end
end

remark = [];
for i = 1:m
    tmp=[];
    for j = 1:N
        tmp=[tmp norm(features(i,:)-u(j,:))];
    end
    [junk, index] = min(tmp);
    remark = [remark; features(i,:) index];
end
end
