```
% Butterworth Filter 
load data
% check if f_s is provided, if not set it below

% Set filter characteristics in digital domain: 
f_c = 125; % cut off freq [Hz]
w_c = 2*pi*f_c/f_s;  % normalized cut off freq [rad/s]
wn = w_c/pi;  % fraction of normalised Nyquist freq.
% omega_c = 2*f_s*tan(w_c/2);  % analog angular cutoff freq [Hz] % Prewarp frequencies
L = 8;   % length of filter = order of filter +1 

% 1-step butterworth filter calculation:
[B, A] = butter(L-1,wn);

% Filter the provided signal:
signalNAME_filtered = filter(B,A,signalNAME);

%% Frequency response of FILTER
% Normalised Nyquist frequency: 
w = pi*(0:floor(N/2))/(N/2); % normalised angular frequency 
f_Hz = (f_s/2)*(0:floor(N/2))/(N/2); % corresponding positive frequency vector (Hz)
  
% Calculating response from the acutal filter: 
H_hat = (B*exp(-1i*(0:L-1)'*w))./ (A*exp(-1i*(0:L-1)'*w));  % corresponds to H_hat(w) = B(e^jnw)/A(e^jnw)

%% Frequency response of SIGNAL <ADJUST AS REQUIRED>
ecg_hat = fft(ecg);  % freq response of given signal 
ecg_filtered_hat  = fft(signalNAME_filtered);  % freq response filtered signal


%% Ploting Frequency Reponse filter (Amplitude & Phase)
figure(1)
clf
subplot(2,1,1) % Plots magnitude 
hold on
plot(f_Hz,mag2db(abs(H_hat)))
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
title('Filter Magnitude freq response')
subplot(2,1,2) % Plots phase  
hold on
plot(f_Hz,angle(H_hat))
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
title('Filter Phase freq response ')

%% Plotting Filtered SIGNAL Output <Adjust as NEEEDED>
figure(2)
clf
% Signal input (plots signal + noise)
subplot(2,1,1)
plot(t,ecg_noise,'b')
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t(end)]);
title('ECG with HF  noise')
hold on
plot(t,ecg,'m')
legend('ecg+noise','ecg')

% ECG filtered output (plots filtered signal + noise)
subplot(2,1,2)
plot(t,ecg_filtered,'g')
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t(end)]);
title('filtered ECG')
hold on
plot(t,ecg,'m')
legend('ecg filtered','ecg')

%% Ploting Frequency Reponse of SIGNAL
figure(3) 
clf
subplot(2,1,1)
hold on
plot(f_Hz,mag2db(abs(ecg_noise_hat(1:floor(N/2)+1))),'b')
%pause
plot(f_Hz,mag2db(abs(ecg_filtered_hat(1:floor(N/2)+1))),'g')
title('Frequency Response of ECG signals')
legend('ecg+noise','ecg filtered')
xlim([0 250])


subplot(2,1,2)
hold on
%pause
plot(f_Hz,mag2db(abs(ecg_hat(1:floor(N/2)+1))),'m')
%pause
plot(f_Hz,mag2db(abs(xHF_hat(1:floor(N/2)+1))),'k')
xlabel('Freq. [Hz]')
ylabel('Amplitude [dB]')
legend('ecg','noise')
xlim([0 250])

% if mag2db doesn't work, use 20*log10
```
