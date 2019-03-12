% Produce the syntax microseismic signal based on Ricker wavelet
function ns = signal()
%------------------------Signal------------------------------------------
%clean signal
s = ricklet();
%subplot 221;
%plot(s);
%title("clean signal");
% noise
n = wgn(1,300,0.25);
n = n-mean(s(:));
n = n/max(abs(n(:)));
%n = normrnd(0,0.4,1,300);
%subplot 222;
%plot(n);
%title("noise");
% noisy signal
ns = awgn(s,1.64,'measured');
%subplot 223;
%plot(ns);
%title("noisy signal");
%grid;  
end