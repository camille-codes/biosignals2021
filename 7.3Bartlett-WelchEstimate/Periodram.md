'''
clear all
N=100; % length of signal 
t=0:N-1;  % time
fs=1;     % normalised sampling frequency
f=2*pi*t/N*fs; % normaised frequency

x=randn(1,N);
X=fft(x);
S=abs(X).^2/N; 

figure(1)
plot(f,S);
xlabel('Normalised Frequency')
ylabel('PSD [a.u.]')
title('Periodogram of WGN')
var(S)
'''
