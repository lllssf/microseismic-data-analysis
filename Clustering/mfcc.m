function ccc = mfcc(s)
[x,fs,bits] = waveread(s);
bank = melbankm(24,256,fs,0,0.5,'m');