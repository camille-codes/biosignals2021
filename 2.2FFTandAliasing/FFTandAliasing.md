```
% ClassEx2_2: square wave aliasing effect
clear all
close

f_s = 1000; %sampling frequency [Hz]
T =1/f_s;   % sampling period [s]
t_end =1.024;  % final time
t = 0:T:t_end;  % time vector (s)
N = length(t); 

f_Hz = (f_s/2)*(0:floor(N/2))/(N/2); % corresponding positive frequency vector (Hz)

f_0 = 30; % [Hz} fundamental frequency of squarewave
x = square(2*pi*f_0*t,50); % squarewave with 50% duty cycle - must have Signal Processing Toolbox installed
x(1)=0;
X = abs(fft(x));  % freq spectrum

%plot signal in time domain
subplot(2,1,1)
plot(t,x)
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
xlim([0 t_end]);
ylim(1.1*[-1 1])

%plot signal in frequency domain
subplot(2,1,2)
plot(f_Hz,X(1:floor(N/2)+1),'b')
xlabel('Freq. [Hz]')
ylabel('Amplitude [a.u.]')



```
