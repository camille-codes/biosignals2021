```
clear all
close all

f_s = 1000; %sampling frequency [Hz]
T =1/f_s;   % sampling period [s]
t_end =10.024;  % final time
t = 0:T:t_end;  % time vector (s)
t_corr = -t_end:T:t_end;  % time vector for autocorrelation (s)
N = length(t); 
f_Hz = (f_s/2)*(0:floor(N/2))/(N/2); % corresponding positive frequency vector (Hz)

% Gaussian white noise
x1 = randn(1,N);
phi_x1 = xcorr(x1,'biased')/N; % autocorrelation
PSD_x1 = 10*(log10(abs(fft(x1)).^2/N));  % PSD [dB]


%Poisson Process
lambda = 10; % Poisson process rate parameter [events/s]
x2 = double(rand(1,N)< lambda/f_s);
phi_x2 = xcorr(x2,'biased')/N; % autocorrelation
PSD_x2 = 10*(log10(abs(fft(x2)).^2/N));  % PSD [dB]

% figure
% plot(10*(log10(abs(fft(phi_x3)))));


%Pink Noise
alpha = 1/0.1; % decay constant for exponential impulse function
h = exp(-alpha*t);  % impulse function
h= h(1:1000);
x3 = conv(x1,h,'same'); 
phi_x3 = xcorr(x3,'biased')/N; % autocorrelation
PSD_x3 = 10*(log10(abs(fft(x3)).^2/N));  % PSD [dB]

%% Plot Signals 
figure(1)
% Gaussian white noise
subplot(3,1,1)
plot(t,x1)
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t_end]);
ylim(2*[-1 1])
title('WGN')

%Poisson Process
subplot(3,1,2)
plot(t,x2)
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t_end]);
ylim(2*[-1 1])
title('Poisson Process')


%Random Walk
subplot(3,1,3)
plot(t,x3)
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t_end]);
%ylim(2*[-1 1])
title('Pink Noise')

%% plot ACFs
figure(2)
subplot(3,1,1)
plot(t_corr,phi_x1)
xlim([-t_end t_end]);
xlabel('time [s]')
ylabel('ACF [a.u.]')
title('WGN')

subplot(3,1,2)
plot(t_corr,phi_x2)
xlim([-t_end t_end]);
xlabel('time [s]')
ylabel('ACF [a.u.]')
title('Poisson Process')

subplot(3,1,3)
plot(t_corr,phi_x3)
xlim([-t_end t_end]);
xlabel('time [s]')
ylabel('ACF [a.u.]')
title('Pink Noise')

%% Plot PSD
figure(3)
subplot(3,1,1)
plot(log2(f_Hz(2:end)),PSD_x1(2:floor(N/2)+1))
xlabel('Freq. [octave]')
ylabel('PSD [dB]')
title('WGN')


subplot(3,1,2)
plot(log2(f_Hz(2:end)),PSD_x2(2:floor(N/2)+1))
xlabel('Freq. [octave]')
ylabel('PSD [dB]')

title('Poisson Process')


subplot(3,1,3)
plot(log2(f_Hz(2:end)),PSD_x3(2:floor(N/2)+1))
xlabel('Freq. [octave]')
ylabel('PSD [dB]')
title('Pink Noise')

```
