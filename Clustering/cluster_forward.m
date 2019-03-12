clear;
clc;

% Initializing
arrivetime = zeros(54,1);
remark = zeros(54,3000);
features = zeros(3000,15);
y = zeros(1,3000);

% get signal
load('doubleforward.mat');
dt=2.0*1.0e-4;
time=0:dt:(length(ob_data_Vz)-1)*dt; 
r_num=54;

datax1=squeeze(ob_data_Vz(1,:,:));
datax11=datax1';
s = datax1;

% 对每个道集的信号进行初至拾取
for i = 1:54
%n = wgn(1,3000,0.25);
%n = n-mean(s(:));
%n = n/max(abs(n(:)));
%n = n/300;
%y = s(i,:)+n;

%y = awgn(s(i,:),3.64,'measured');
%datax11(:,i) = y';

y = s;

yw_max = zeros(1,3000);
y_w = zeros(1,21);
for  t=1:3000
    if t<=10 
        for j = 1:21
        y_w(j) = y(t+j);
        end
    elseif t>=2990
        for k = 1:21
        y_w(k) = y(t-k);
        end
    else 
        for l = 1:21
        y_w(l) = y(t-10+l);
        end
    end
    yw_max(t) = max(y_w);
    yw_min(t) = min(y_w);
    yw_M(t) = mean(y_w);                    %平均值
    yw_peak(t) = yw_max(t) - yw_min(t);     %峰-峰值
    yw_av(t) = mean(abs(y_w));              %绝对值的平均值
    yw_var(t) = var(y_w);                   %方差
    yw_std(t) = std(y_w);                   %标准差
    yw_kurt(t) = kurtosis(y_w);             %峭度
    yw_skew(t) = skewness(y_w);             %偏度
    yw_rms(t) = rms(y_w);                   %均方根
    yw_S(t) = yw_rms(t)/yw_av(t);           %波形因子
    yw_C(t) = yw_peak(t)/yw_rms(t);         %峰值因子
    yw_I(t) = yw_peak(t)/yw_av(t);          %脉冲因子
    yw_xr(t) = mean(sqrt(abs(y_w)))^2;
    yw_L(t) = yw_peak(t)/yw_xr(t);          %裕度因子
    yw_E(t) = sum(y_w.^2);                  %能量
end
%yw_skew(isnan(yw_skew)) = 0;

features(:,1) = yw_max'; %495  433
%features(:,2) = yw_min'; %588/486  417/407
features(:,3) = yw_M; %580/491  425
%features(:,4) = yw_peak; %546  408
%features(:,5) = yw_av; %579  415
%features(:,6) = yw_var; %560  467
%features(:,7) = yw_std; %546  410
%features(:,8) = yw_kurt; %NAN
%features(:,9) = yw_skew; %NAN
%features(:,10) = yw_rms; %578  413
%features(:,11) = yw_S; %NaN
%features(:,12) = yw_C; %NaN
%features(:,13) = yw_I; %NaN
%features(:,14) = yw_L; %NaN
%features(:,15) = yw_E; %587  456
remark(i,:) = kmeans(features,3);
for t1 = 1:3000
    if (remark(i,t1+1)-remark(i,t1))~=0
        break;
    end
end

arrivetime(i,1) = t1;
    
end
arrivetime = arrivetime'.*dt;

set(gcf,'outerposition',get(0,'screensize')); 
wigb(datax11,1,1:54,time)
	set(gca,'FontName','Times New Roman','FontSize',18);
	title(['Original Seismic Record'],'FontName','Times New Roman','Fontsize',18,'Fontweight','Bold');
	xlabel('{{\it{Trace}}}','FontName','Times New Roman','FontSize',18);
	ylabel('{{\it{Time}}/s}','FontName','Times New Roman','FontSize',18);
hold on;
scatter(1:54,arrivetime(1:54),'ro')
%subplot 122;
%plot(s(54,:));
%hold on;
%tt =arrivetime(54,1);
%plot(tt,s(54,tt),'ro');


