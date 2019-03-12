% Produce the syntax microseismic signal based on Ricker wavelet
function s = ricklet()
% 模拟信号，雷克子波
f = 30;
fs = 1024;
s1 = zeros(1,300);
s2 = zeros(1,300);
s3 = zeros(1,300);
s = zeros(1,300);
for t = 1:300
  s1(t) = (1-2*(pi*70*(t-100)/fs).^2) .*exp(-(pi*70*(t-100)/fs).^2);
  s2(t) = (1-2*(pi*70*(t-150)/fs).^2) .*exp(-(pi*70*(t-150)/fs).^2);
  s3(t) = (1-2*(pi*70*(t-200)/fs).^2) .*exp(-(pi*70*(t-200)/fs).^2);
  s(t) = s1(t) + (2^0.5)*s2(t) + (3^0.5)*s3(t);
end
%plot(s);
end