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
hC = exp(-(t-t(N/2)).^2/2/sigma^2).*cos(2*pi*F0*t);  hC = hC/norm(hC); % filter in cos phase
hR = randn(1,N);  % random filter - for comparison 
hR=hR-(hR*hS')*hS; hR = hR/norm(hR);


% stimulus generation
tL = (1:R)*T; % time for signal
x=randn(R,1);  % Gaussian white noise stimulus 
%x = 100*[babble(1:2:end) babble(2:2:end)]';

% neural response
% step 1: filtered stimulus
vS = filter(hS,1,x)/100; 
vC = filter(hC,1,x)/100; 
vR = filter(hR,1,x)/100; 


% step 2: response rate
A = 600; % approx maximum response rate  [spike/s]
Th = 1e-3; % threshold for spike
SF = @(xx) A*max(xx-Th,0).^2*1000; %one-sided spiking non-liearity 
%SF = @(xx) A*max(abs(xx)-Th,0).^2*1000;% quadratic spiking non-linearity
y = SF(vS);  % response rate

% create vectors for plotting non-linearity
vplot = linspace(-0.05,0.05);
yplot = SF(vplot); % 

% step 3: conversion to spike train Poisson process
spk = (y*T > rand(R,1));

% plot model
figure(1), clf
subplot(1,2,1)
hold on
plot(t*1000,hS,'LineWidth',2)
xlabel('time [ms]')
ylabel('amplitude [a.u]')
title('system filter')

subplot(1,2,2)
hold on
plot(1000*vplot,yplot,'LineWidth',2)
xlabel('generator signal [mV]')
ylabel('response [spikes/s]')
title('system nonlinearity')

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

figure(3), clf
v_ref = vR;
hold on
plot(vS,v_ref,'k.')
plot(vS(spk),v_ref(spk),'m.')
axis equal
xlabel('Generator Signal with true filter hS')
ylabel('Generator Signal with orthogonal filter h_ref')


figure(4), clf
hold on
NN =length(vR);
scatter3(vR(1:100:NN),vC(1:100:NN),vS(1:100:NN),'k.')
scatter3(vR(spk),vC(spk),vS(spk),'m.')
axis equal
view(2)

xlabel('Gen. Sig. with Random filter ')
ylabel('Gen. Sig. with Cosine filter')
zlabel('Gen. Sig. with System filter')

figure(5), clf
hold on
subplot(3,1,1)
plot(t*1000,hS,'LineWidth',2)
title('system  filter - sine phase')

subplot(3,1,2)
plot(t*1000,hC,'LineWidth',2)
subplot(3,1,3)
plot(t*1000,hR,'LineWidth',2)
title('orthogonal reference filter - random')
xlabel('time [ms]')
ylabel('amplitude [a.u]')
title('orthogonal reference filter - cosine phase')


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



%% Nonlinearity estimation
B = 30; % number of bins in histogram for nonlinearity estimation
v_STA = filter(STA,1,x)/100;  % project data onto STA = prediction of generator signal
v_STA_spk = v_STA(spk);  % find those generator signals giving a spike

[Nfull,edges] = histcounts(1000*v_STA,B); % histogram for full stimulus distribution
[Nspike,~] = histcounts(1000*v_STA_spk,edges); % histogram for spiking triggered distribution

SF_STA = Nspike./Nfull/T;

centres =(edges(1:end-1)+edges(2:end))/2; % ca


figure(6)
subplot(3,1,1)
bar(centres,Nfull,'k')
xlabel('generator signal [mV]')
ylabel('occurence (#)')
title('Full stimulus distribution of generator signal')


subplot(3,1,2)
bar(centres,Nspike,'m')
xlabel('generator signal [mV]')
ylabel('occurence (#)')
title('Spike triggered stimulus distribution of generator signal')


subplot(3,1,3)
plot(centres, SF_STA,'or-')
xlabel('generator signal [mV]')
ylabel('response [spikes/s]')
title('Estimated system nonlinearity')

figure(1)
subplot(1,2,2)
hold on 
plot(centres, SF_STA,'or-')

```
