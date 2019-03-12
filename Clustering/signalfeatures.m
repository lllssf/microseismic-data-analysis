function [y,features] = signalfeatures(s)


% --------------------Window length = 21
y = s;

%---------------------Initializing features
yw_max = zeros(1,300);
yw_min = zeros(1,300);
yw_M = zeros(1,300);                   
yw_peak = zeros(1,300);     
yw_av = zeros(1,300);              
yw_var= zeros(1,300);                
yw_std = zeros(1,300);                 
yw_kurt = zeros(1,300);            
yw_skew = zeros(1,300);           
yw_rms = zeros(1,300);                  
yw_S = zeros(1,300);         
yw_C = zeros(1,300);        
yw_I = zeros(1,300);          
yw_xr = zeros(1,300);
yw_L = zeros(1,300);        
y_E = zeros(1,300);

%-------每一帧的feature
for  t=1:300
    if t<=10 
        for i = 1:21
        y_w(i) = y(t+i);
        end
    elseif t>=290
        for i = 1:21
        y_w(i) = y(t-i);
        end
    else 
        for i = 1:21
        y_w(i) = y(t-10+i);
        end
    end
    if t<=0 || t>=900
    yw_max(t) = 0;
    yw_min(t) = 0;
    yw_M(t) = 0;                    %平均值
    yw_peak(t) = 0;     %峰-峰值
    yw_av(t) = 0;              %绝对值的平均值
    yw_var(t) = 0;                   %方差
    yw_std(t) = 0;                   %标准差
    yw_kurt(t) = 0;             %峭度
    yw_skew(t) = 0;             %偏度
    yw_rms(t) =0;                   %均方根
    yw_S(t) = 0;           %波形因子
    yw_C(t) = 0;         %峰值因子
    yw_I(t) = 0;          %脉冲因子
    yw_xr(t) = 0;
    yw_L(t) = 0;          %裕度因子
    yw_E(t) = 0;
    
    else 
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
        
end
%subplot 224
%plot(yw_M)
%axis([0 250 -1 1]);

features = zeros(300,5);
features(:,1) = yw_max';
%features(:,2) = yw_min';
%features(:,2) = yw_M;
%features(:,3) = yw_peak;
%features(:,5) = yw_av;
%features(:,2) = yw_var;
%features(:,7) = yw_std;
%features(:,8) = yw_kurt;
%features(:,9) = yw_skew;
%features(:,5) = yw_rms;
%features(:,3) = yw_S;
%features(:,4) = yw_C;
%features(:,5) = yw_I;
%features(:,14) = yw_L;
%features(:,15) = yw_E;



end