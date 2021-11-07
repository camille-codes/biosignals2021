```
clear all


%% Simulation of system to gather data
N=50; % dimensionality for filter
R= 10^5;  % length of stimulus sequence
fs=500;  % sampling rate [Hz]
T=1/fs;   % sample time step
t=(1:N)*T;  % time for filter [s]

% filter specification
sigma = 0.020;  % duration of Gaussian filter envelope [s]
F0 = 50;  % preferred frequency of filter
hS = exp(-(t-t(N/2)).^2/2/sigma^2).*sin(2*pi*F0*t); hS = hS/norm(hS);% filter in sin phase
hC = exp(-(t-t(N/2)).^2/2/sigma^2).*cos(2*pi*F0*t);  hC = hC/norm(hC); % filter in cos phase
hR = randn(1,N);  % random filter - for comparison 
hR=hR-(hR*hS')*hS -(hR*hC')*hC; hR = hR/norm(hR);
hR2 = randn(1,N);  % 2nd random filter - for comparison 
hR2=hR2-(hR2*hS')*hS -(hR2*hC')*hC - (hR2*hR')*hR; hR2 = hR2/norm(hR2);



% stimulus generation
tL = (1:R)*T; % time for signal
x=randn(R,1);  % Gaussian white noise stimulus 
%x = 100*[babble(1:2:end) babble(2:2:end)]';

% neural response
% step 1: filtered stimulus
vS = filter(hS,1,x)/100; 
vC = filter(hC,1,x)/100; 
vR = filter(hR,1,x)/100; 
vR2 = filter(hR2,1,x)/100; 


% step 2: response rate
A = 300; % approx maximum response rate  [spike/s]
Th = 0; % threshold for spike
%SF = @(xx) A*max(xx-Th,0).^2*1000; %one-sided spiking non-liearity 
SF = @(xx) A*max(abs(xx)-Th,0).^2*1000;% quadratic spiking non-linearity
y = 10*sqrt(SF(vS)+SF(vC));  % response rate

% create vectors for plotting non-linearity
vplot = linspace(-0.05,0.05,101);
yplot = 10*sqrt(SF(vplot)+SF(vplot')); % 
yplot1 = yplot(:,51)';
yplot2 = yplot(51,:);

% step 3: conversion to spike train Poisson process
spk = (y*T > rand(R,1));

% plot model
figure(1), clf
subplot(2,3,1)
hold on
plot(t*1000,hS,'b','LineWidth',2)
xlabel('time [ms]')
ylabel('amplitude [a.u]')
title('system filter - sine phase')

subplot(2,3,2)
hold on
plot(1000*vplot,yplot1,'b','LineWidth',2)
xlabel('vS')
ylabel('response [spikes/s]')
title('1D nonlinearity - sine phase')


subplot(2,3,4)
hold on
plot(t*1000,hC,'g','LineWidth',2)
xlabel('time [ms]')
ylabel('amplitude [a.u]')
title('system filter - cosine phase')

subplot(2,3,5)
hold on
plot(1000*vplot,yplot2,'g','LineWidth',2)
xlabel('vC')
ylabel('response [spikes/s]')
title('1D nonlinearity - cosine phase')

subplot(2,3,3)
hold on
contour(1000*vplot,1000*vplot,yplot)
plot(0*vplot,1000*vplot,'b--')
plot(1000*vplot,0*vplot,'g--')
xlabel('vS ')
ylabel('vC')
cc= colorbar;
cc.Label.String = 'spike rate [Hz]';
axis square
title('2D Spiking Nonlinearity') 

% plot model signals 
figure(2) , clf 
subplot(4,1,1)
plot(tL*1000,x,'LineWidth',2)
xlabel('time [ms])')
ylabel('amplitude [a.u.]')
title('stimulus signal')
tmax = 500;
xlim([0 tmax])

subplot(4,1,2)
hold on
plot(tL*1000,1000*vS,'b','LineWidth',2)
xlabel('time [ms])')
ylabel('vS [mV]')
title('output sine phase filter')
xlim([0 tmax])

subplot(4,1,3)
hold on
plot(tL*1000,1000*vC,'g','LineWidth',2)
xlabel('time [ms])')
ylabel('vC [mV]')
title('output cosine phase filter')
xlim([0 tmax])

subplot(4,1,4)
hold on
plot(tL*1000,y,'LineWidth',2)
plot(tL(spk)*1000,y(spk),'mo')
xlabel('time [ms])')
ylabel('response [spike/s]')
title('neural response')
xlim([0 tmax])

figure(3), clf
subplot(1,2,1)
hold on
plot(vS,vC,'k.')
plot(vS(spk),vC(spk),'m.')
axis equal
xlabel('vS [a.u]')
ylabel('vC [a.u]')
title('Stimulus space projected onto sine and cosine filters')
subplot(1,2,2)
hold on
plot(vR,vR2,'k.')
plot(vR(spk),vR2(spk),'m.')
axis equal
xlabel('vR1 [a.u]')
ylabel('vR2 [a.u]')
title('Stimulus space projected onto sine and cosine filters')


figure(4), clf
hold on
sb =5;
scatter3(vR(1:sb:end),vC(1:sb:end),vS(1:sb:end),'k.')
scatter3(vR(spk),vC(spk),vS(spk),'m.')
axis equal
view(2)

xlabel('vR [a.u]')
ylabel('vC [a.u]')
zlabel('vS [a.u]')
title('Stimulus space projected onto sine, cosine and random filters')


return

%% Spike triggered covariance
%spk = double(spk);

for n =1:N
    X(n,:) = x(N+1-n:end-n+1)';
    XS(n,:) = (double(spk(N:end)).* x(N+1-n:end-n+1))';
end
spikecount= sum(spk(N:R)); 
CS=XS*XS'/(spikecount);
C = X*X'/(R-N);

DC = CS-C;

[V D] = eigs(DC,N); 

figure(5)
plot(diag(D),'k.')
xlabel('eigenvalue index')
ylabel('eigenvalue')
title('STC eigenvalues')

STC1 = flipud(V(:,1));
STC2 = flipud(V(:,2));

figure(1)
subplot(2,3,1)
plot(t*1000,STC1,'r--','LineWidth',2)
xlabel('time [ms]')
ylabel('amplitude [a.u]')
legend('true','STA')

subplot(2,3,4)
plot(t*1000,STC2,'r--','LineWidth',2)
xlabel('time [ms]')
ylabel('amplitude [a.u]')
legend('true','STA')

%return
%% Nonlinearity estimation
B = 30; % number of bins in histogram for nonlinearity estimation
v_STC1 = filter(STC1,1,x)/100;  % project data onto STC1 = prediction of generator signal
v_STC2 = filter(STC2,1,x)/100;  % project data onto STC2 = prediction of generator signal
v_STC1_spk = v_STC1(spk);  % find those generator signals giving a spike
v_STC2_spk = v_STC2(spk);  % find those generator signals giving a spike

[Nfull,edges1,edges2] = histcounts2(1000*v_STC1,1000*v_STC2,B); % histogram for full stimulus distribution
[Nspike,~] = histcounts2(1000*v_STC1_spk,1000*v_STC2_spk,edges1,edges2); % histogram for spiking triggered distribution

SF_STC = Nspike./Nfull/T; % nonlinearity estimate.

centres1 =(edges1(1:end-1)+edges1(2:end))/2; % 
centres2 =(edges2(1:end-1)+edges2(2:end))/2; % 


%plot 2D nonlinearity
figure(1)
subplot(2,3,6)
imagesc(SF_STC)
axis square 
xlabel('v1')
ylabel('v2')
colorbar

% plot 1D effective nonlinearities
figure(1)
subplot(2,3,2)
hold on 
plot(centres1, SF_STC(B/2+1,:),'or-')

subplot(2,3,5)
hold on 
plot(centres2, SF_STC(:, B/2+1),'or-')

```
