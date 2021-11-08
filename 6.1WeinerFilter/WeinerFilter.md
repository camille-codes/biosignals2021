clear all
%close all

load ECGdata1000

%f_s = 1000; %sampling frequency [Hz] loaded from ECGdata1000
T =1/f_s;   % sampling period [s]
N = length(t); 
A_HF = 0.04; % High Frequency Noise Amplitude
A_LF = 2; % Low Frequency Noise Amplitude
A_50 = 0.1; % % 50 Hz power line noise Amplitude



% ECG signal
xHF = A_HF*max(ecg)*randn(1,N); % High Frequency Gaussian white noise
xLF = A_LF*max(ecg)*t/t(end); % Low Frequency noise
xLF = xLF -mean(xLF); 
x50 =  A_50*sin(2*pi*50*t); % 50 Hz noise

% choose noise type
ns =1; 
if ns ==1
    x = xHF;  
elseif ns ==2
    x = xLF;  
elseif ns == 3
    x = x50;  
else
    return
end
    
    
        
ecg_noise = ecg+x; % ECG signal with high and low freq noise


%% Characterize PSD of ideal signal and noise
%normalised Nyquist frequency 
M = round(f_s/2/li)+1; % length of one cycle of ECG
f_Hz_1 = (f_s/2)*(0:floor(M/2))/(M/2); % corresponding positive frequency vector (Hz)



ecg_hat =fft(ecg(1:M));  % freq response ecg
x_hat = fft(x(1:M));  % freq response  noise
x_hat_id = x_hat; 
if ns==1
    x_hat_id = sqrt(M)*A_HF*max(ecg)*ones(1,M); % freq response of HF Gaussian white noise
end

% plot PSDs
figure(4) 
clf

hold on
plot(f_Hz_1,10*log10(abs(ecg_hat(1:floor(M/2)+1)).^2/N),'m')
plot(f_Hz_1,10*log10(abs(x_hat(1:floor(M/2)+1)).^2/N),'k')
plot(f_Hz_1, 10*log10(abs( x_hat_id(1:floor(M/2)+1)).^2/N),'r')
xlabel('Freq. [Hz]')
ylabel('Amplitude [dB]')
title ('PSD of ecg signal and noise')

legend('ecg','noise empirical','noise theoretical')

xlim([0 250])



%% Specify filter impulse response

L= ;   % length of filter

A = [1];  % a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%B = % input inpulse response function in terms of n and w_c.
 % - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
 
Ww =  ;  % fft of the Wiener filter weights
ww = ifftshift(ifft(Ww)); % Full length Wiener filter, shifted to be symmetric about centre
centre = ceil(M/2)+(-floor(L/2):floor(L/2));  % central portion of length L to extract
B = ww(centre);  % final filter weights
% plot filter
figure(5)
clf
hold on
plot(ww,'b')
plot(centre,B,'g')
title('Wiener filter')
xlabel('time [samples]')
ylabel('amplitude [a.u.]')
legend('full filter', 'central filter L samples')



%% Filter signal 
ecg_filtered = filter(B,A,ecg_noise);


%% Frequency response of filter
%normalised Nyquist frequency 
w = pi*(0:floor(N/2))/(N/2); % normalised angular frequency 
f_Hz = (f_s/2)*(0:floor(N/2))/(N/2); % corresponding positive frequency vector (Hz)
w = w + 1e-6; % too avoid division by zero  

% calculate response from acutal filter
H_hat = B*exp(-1j*(0:L-1)'*w);  % corresponds to H_hat(w) = sum h[n] e^(-jwn)




%% Frequency response of signal 
ecg_hat = fft(ecg);  % freq response ecg
x_hat = fft(x);  % freq response HF noise
ecg_noise_hat = fft(ecg_noise);  % freq response input signa
ecg_filtered_hat  = fft(ecg_filtered);  % freq response filtered signal

%% Ploting Frequency Reponse filter

figure(6)
clf
hold on
plot(f_Hz,20*log10(abs(H_hat)),'k')
plot(f_Hz_1,20*log10(abs(Ww(1:length(f_Hz_1)))),'r--')

xlabel('Freq. [Hz]')
ylabel('Amplitude [dB]')
title('Filter freq response magnitude')
legend('true','ideal')
 



%% Plotting Filtered Output
figure(7)
clf
% ECG input 
subplot(2,1,1)
plot(t,ecg_noise,'b')
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t(end)]);
title('ECG with HF  noise')
hold on
plot(t,ecg,'m')
legend('ecg+noise','ecg')
xlim([0 10])


% ECG filtered output
subplot(2,1,2)
plot(t,ecg_filtered,'g')
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t(end)]);
title('filtered ECG')
hold on
plot(t,ecg,'m')
legend('ecg filtered','ecg')
xlim([0 10])

%% Ploting Frequency Reponse signal
figure(8) 
clf
subplot(2,1,1)
hold on
plot(f_Hz,20*log10(abs(ecg_noise_hat(1:floor(N/2)+1))),'b')
%pause
plot(f_Hz,20*log10(abs(ecg_filtered_hat(1:floor(N/2)+1))),'g')
title('Frequency Response of ECG signals')
legend('ecg+noise','ecg filtered')
xlim([0 250])


subplot(2,1,2)
hold on
%pause
plot(f_Hz,20*log10(abs(ecg_hat(1:floor(N/2)+1))),'m')
%pause
plot(f_Hz,20*log10(abs(x_hat(1:floor(N/2)+1))),'k')
xlabel('Freq. [Hz]')
ylabel('Amplitude [dB]')


legend('ecg','noise')
xlim([0 250])


 

