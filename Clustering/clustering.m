clear
clc
%function arrivetime = clustering()
% -------------------------clean signal arrive time
%clean signal
[y,features] = signalfeatures(ricklet());
clean_remark = kmeans(features,2);
y_re = y;
arrivetime1 = ones(3,1);
i = 1;
for t2 = 1:299
    if (clean_remark(t2+1)-clean_remark(t2)) ~= 0
        arrivetime1(i) = t2;
        i = i+1;
    end
end
subplot 211
plot(y);
title("clean signal");
hold on;
for i  = 1:2:length(arrivetime1)-1
plot(arrivetime1(i),y(arrivetime1(i)),'ro');
end


%------------------------K-Means
%[y,features,feature_index] = featureselect();
[y,features] = signalfeatures(signal());
remark = kmeans(features,2);
y_re = y;
arrivetime2 = ones(3,1);
i = 1;
for t2 = 1:299
    if (remark(t2+1)-remark(t2)) ~= 0
        arrivetime2(i) = t2;
        i = i+1;
    end
end
subplot 212
plot(y);
title("noisy signal");
hold on;
for i  = 1:2:length(arrivetime2)-1
plot(arrivetime2(i),y(arrivetime2(i)),'ro');
end


%end
