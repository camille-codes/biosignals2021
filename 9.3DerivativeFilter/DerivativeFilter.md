```
clear all
%close all

load ECGdata360
val=  ecg1;
%val=  ecg2;
%val=  ecg3;
%val=  ecg4;

ecg = val/max(val);
f_s = 360; %sampling frequency [Hz] loaded from ECGdata1000
T =1/f_s;   % sampling period [s]
t = (1:length(ecg))*T;
A_HF = 0.0; % High Frequency Noise Amplitude
ecg  = ecg + A_HF*randn(1,length(t));



%% Filtering 
x = ecg; 
d = [diff(x) 0]; % derivative filter
d2 = d.^2;  % sqaure of derivative
N = 12;  % length of filter for g1 sum
h1 = fliplr(0:N); % filter for g1 sum 
g1 = filter(h1,1,d2);

M = 12; % length of filter for moving average to get g
h2 = ones(1,M)/M;  % filter for moving average to get g
g = filter(h2,1,g1);

gth = 0.04; % threshold for  QRS complex detection

suprath = (g> gth); % suprathreshold point
L = 6; % window for finding peaks
icnt = 0;
for l =1:length(suprath)
    start = max(l-L,1);
    stop = min(length(suprath),l+L);
    if suprath(l) & (g(l) == max(g(start:stop)))
        icnt = icnt +1;
        peakss(icnt) = l;
    end
end
      

        



%% Plotting Filtered Output
figure(1), clf
% ECG input 
subplot(5,1,1)
plot(t,x)
%xlabel('Time [s]')
ylabel('x')
xlim([0 t(end)]);
title('x= ECG')
hold on

subplot(5,1,2)
plot(t,d)
%xlabel('Time [s]')
ylabel('d')
xlim([0 t(end)]);
title('d')
hold on

subplot(5,1,3)
plot(t,d2)
%xlabel('Time [s]')
ylabel('d^2')
xlim([0 t(end)]);
title('d^2')
hold on

subplot(5,1,4)
plot(t,g1)
%xlabel('Time [s]')
ylabel('g1')
xlim([0 t(end)]);
title('g1')
hold on

subplot(5,1,5)
plot(t,g)
%xlabel('Time [s]')
ylabel('g')
xlim([0 t(end)]);
title('g')
hold on
plot(t(suprath),g(suprath),'g.')
plot(t(peakss),g(peakss),'r*')



delay = (N+M)/2;
peakss = peakss -delay;
suprath = suprath(delay:end);

subplot(5,1,1)
hold on
plot(t(suprath),x(suprath),'g.')
plot(t(peakss),x(peakss),'r*')
legend('signal','QRS complex','peaks')




```
