```
clear all
%% Load vowel data
% vIY is a waveform for vowel 'IY'
% vIY_AH is a waveform for 'IY' and then 'AH'
% vall is a waveform for all 10 vowels, one after the other
% str.vowels contains name of the vowels

load 'VowelData2.mat'

N = length(vIY);

%% Obtain matched filter for 'IY' 
figure(1), clf
subplot(2,1,1)
hold on
plot(vIY)
xlabel('Time [samples]')
ylabel('Amplitude [a.u]')
title('IY - time in samples ')
%s =1:100;
s = 99:199;  % choose appropriate range for samples for matched filter
h = vIY(s); % wavefrom for extract first pitch period
h = fliplr(h'); % reverse the waveform in time to get the matched filter
plot(s,vIY(s),'r')
xlim([0 5000])
subplot(2,1,2)
plot(t,vIY)
xlabel('Time [s]')
ylabel('Amplitude [a.u]')
title('IY- full signal ')

return
%% Filter 'vIY' signal and plot
A_HF = 0.05; % amplitude of high frequency noise
x = vIY +A_HF*randn(N,1);
y = conv(x,h,'same');
sound(x,fs)

figure(2), clf
subplot(2,1,1)
plot(t, x)
xlabel('t [s]')
ylabel('Amplitude [a.u]')
title('IY with white noise')
subplot(2,1,2)
plot(t, y)
xlabel('t [s]')
ylabel('Amplitude [a.u]')
title('matched IY-filter output ')

return
%% Filter 'vIY+AH' signal and plot
A_HF = 0.05; % amplitude of high frequency noise
x = vIY_AH +A_HF*randn(N,1);
y = conv(x,h,'same');
sound(x,fs)

figure(3), clf
subplot(2,1,1)
plot(t, x)
xlabel('t [s]')
ylabel('Amplitude [a.u]')
title('IY_AH with white noise')
subplot(2,1,2)
plot(t, y)
xlabel('t [s]')
ylabel('Amplitude [a.u]')
title('matched IY-filter output')
return 
%% Filter 'vall' signal and plot
A_HF = 0.05; % amplitude of high frequency noise
x = vall +A_HF*randn(N,1);
y = conv(x,h,'same');
sound(x,fs)

figure(4), clf
subplot(2,1,1)
plot(t, x)
xlabel('t [s]')
ylabel('Amplitude [a.u]')
title('all vowel with white noise')
subplot(2,1,2)
plot(t, y)
xlabel('t [s]')
ylabel('Amplitude [a.u]')
title('matched IY-filter output')

```
