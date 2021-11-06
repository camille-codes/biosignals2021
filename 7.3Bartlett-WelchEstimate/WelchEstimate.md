```
clear all
close all

file = 'C:\Users\hmeffin\Home\Teaching\2021\BMEN90035 Biosignal Processing\Lectures\Vowels\Data\RH\HAD_00.wav';
[y, fs] = audioread(file);

T =1/fs;   % sampling period [s]
N = length(y); 
t=(0:N-1)*T;
it = 0:N-1;

%% plot  signal
figure(1), clf 
subplot(2,1,1)
%plot(it,y)
plot(t,y)
xlabel('Time [s]')
ylabel('Amplitude SPL [a.u.]')
title('Vowel Signal') 

%% Extract vowel
t_start = T;
t_end =t(end);
it_start = round(t_start/T);
it_end= round(t_end/T);
y = y(it_start:it_end);
N = length(y); 
t=(0:N-1)*T;
it = 0:N-1;
hold on
%plot(it_start+it,y)
plot(t_start+t,y)

%% Plot PSD from full signal
f_Hz2 = (fs/2)*(0:floor(N/2))/(N/2); % corresponding positive frequency vector (Hz)
S_long = 2*T/N*abs(fft(y)).^2;  % PSD using whole signal
figure(2), clf
subplot(2,2,1)
plot(f_Hz2,10*log10(S_long(1:floor(N/2)+1)))
fmax  = 5000; % max plotting freq [Hz]
xlim([0 fmax])
xlabel('Freq. [Hz]')
ylabel('PSD [dB/Hz]')
title('PSD full signal') 
ylim([-100 -50])

%% Specify parameters for Welch Estimate
M = 1*430+1; % window length - choose to be odd
w = rectwin(M); % window function to apply directly to signal
w = w/norm(w); % normalise to unit energy
L=floor(N/M); % number of segments

Y = reshape(y(1:L*M),M,L)';  % split signal into L segments of length M
Y = Y.* repmat(w',L,1);  % multiply by window function

phi_sum= zeros(1,2*M-1); % initialize sum for ACF
S_sum =  zeros(1,M); % initialize sum for PSD 
figure(1)
f_Hz = (fs/2)*(0:floor(M/2))/(M/2); % corresponding positive frequency vector (Hz)
Rp = 1;
for ll =1:L
    phi_sum = phi_sum + xcorr(Y(ll,:)); % update ACF
    S_sum = S_sum + 2*T*abs((fft(Y(ll,:)))).^2; % update PSD
        figure(2)
        subplot(2,2,2)
        plot(f_Hz,10*log10(S_sum(1:floor(M/2)+1)/ll)) % plot running estimate of S = S_sum/ll
        xlabel('Freq. [Hz]')
        ylabel('PSD [dB/Hz]')
        title(['PSD l=' num2str(ll)]) 
        ylim([-100 -50])
        xlim([0 fmax])
        figure(1)
        subplot(2,1,2)
        plot((-M+1:M-1)*T,phi_sum/ll)  % plot running estimate of phi = phi_sum/ll
        xlabel('Time. [s]')
        ylabel('ACF [a.u.]')
        title(['ACF l=' num2str(ll)]) 
        %pause
end


%% Checking answer
figure(2)
% check by convolving in frequency domain.
NN=N;
om = pi*(-floor(NN/2):floor(NN/2))/(NN/2); % normalised angular frequency 
om = om + 1e-6; % too avoid division by zero  
f_Hz2 = (fs/2)*(0:floor(NN/2))/(NN/2); % corresponding positive frequency vector (Hz)
%w_hat = 1/MM/sqrt(M)*(sin(om*MM/2)./sin(om/2)).^2; %  freq response for Bartlett window
w_hat = T/M*(sin(om*M/2)./sin(om/2)).^2; %  freq response for Bartlett window

S_long = T/N*abs((fft(y(1:NN)))).^2;  % PSD using whole signal
S_B2 = 2*conv(w_hat,S_long,'same');
subplot(2,2,3)
plot(f_Hz2,10*log10(S_B2(1:floor(NN/2)+1)))
xlim([0 fmax])
xlabel('Freq. [Hz]')
ylabel('PSD [dB/Hz]')
title('Convolutional PSD')
ylim([-100 -50])

%check using Matlab built-in function
subplot(2,2,4)
pwelch(y,w,0,[],fs,'onesided')
xlim([0 fmax/1000])
ylim([-100 -50])

figure(2)
subplot(2,2,1)
hold on
plot(f_Hz,10*log10(S_sum(1:floor(M/2)+1)/ll),'r')

```
