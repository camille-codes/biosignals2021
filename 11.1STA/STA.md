```
clear all
load Babble

%% Simulation of system to gather data
N=500; % dimensionality for filter
R= length(babble);  % length of stimulus sequence
fs=5000;  % sampling rate [Hz]
T=1/fs;   % sample time step
t=(1:N)*T;  % time for filter [s]

% filter specification
sigma = 100/fs;  % duration of Gaussian filter envelope [s]
F0 = 50;  % preferred frequency of filter
hS = exp(-(t-t(N/2)).^2/2/sigma^2).*sin(2*pi*F0*t); hS = hS/norm(hS);% filter in sin phase
%hR=hR-(hR*hS')*hS; hR = hR/norm(hR);

% stimulus generation
tL = (1:R)*T; % time for signal
x=randn(R,1);  % Gaussian white noise stimulus 
%x = 100*[babble(1:2:end) babble(2:2:end)]'; % multispeaker babble stimulus

sound(x(1:5*fs),fs);

% neural response 
% step 1: filtered stimulus
vS = filter(hS,1,x)/100;  

% step 2: response rate
A = 600; % approx maximum response rate  [spike/s]
Th = 1e-3; % threshold for spike
SF = @(xx) A*max(xx-Th,0).^2*1000; %one-sided spiking non-liearity 
SF = @(xx) A*max(abs(xx)-Th,0).^2*1000;% quadratic spiking non-linearity
y = SF(vS);  % response rate

% create vectors for plotting non-linearity
vplot = linspace(-0.03,0.03);
yplot = SF(vplot);  % 


% step 3: conversion to spike train Poisson process
spk = (y*T > rand(R,1));

% plot model
figure(1), clf
subplot(1,2,1)
hold on
plot(t*1000,hS,'LineWidth',2)
xlabel('time [ms]')
ylabel('amplitude [a.u]')
title('ANF filter')

subplot(1,2,2)
hold on
plot(1000*vplot,yplot,'LineWidth',2)
xlabel('generator signal [mV]')
ylabel('response [spikes/s]')
title('ANF nonlinearity')

figure(2) , clf 
subplot(3,1,1)
plot(tL*1000,x,'LineWidth',2)
xlabel('time [ms])')
ylabel('amplitude [a.u.]')
title('stimulus signal')
tmax = 500;
xlim([0 tmax])

subplot(3,1,2)
hold on
plot(tL*1000,1000*vS,'LineWidth',2)
xlabel('time [ms])')
ylabel('generator signal [mV]')
title('generator signal')
xlim([0 tmax])

subplot(3,1,3)
hold on
plot(tL*1000,y,'LineWidth',2)
plot(tL(spk)*1000,y(spk),'mo')
xlabel('time [ms])')
ylabel('response [spike/s]')
title('neural response')
xlim([0 tmax])

return

%% Spike triggered Average
%spk = double(spk);
spikecount= sum(spk(N:R)); 
STA = xcorr(double(spk),x,N)/spikecount;
STA = STA(end-N+1:end);
STA = STA/norm(STA);

figure(1)
subplot(1,2,1)
plot(t*1000,STA,'r--','LineWidth',2)
xlabel('time [ms]')
ylabel('amplitude [a.u]')
legend('true','STA')
```
