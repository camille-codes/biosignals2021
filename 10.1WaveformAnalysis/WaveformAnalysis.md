![Image of Yaktocat](https://github.com/camille-codes/biosignals2021/blob/main/10.1WaveformAnalysis/10.1.png)
```
clear all
%close all

load('ECGdataPVC.mat')

T =1/f_s;   % sampling period [s]
t = (1:numel(data))*T;

%% Plot Data
% ECG input 
figure(1)
plot((1:size(data,2))*T,data(1,:))
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
title('ECG with PVC')
hold on

return

%% Filtering & Plotting
x = reshape(data',numel(data),1)'; 
x = x/max(x);
d = [diff(x) 0]; % derivative filter
d2 = d.^2;  % sqaure of derivative
N = 12;  % length of filter for g1 sum
h1 = 0:N; % filter for g1 sum
%g1 = conv(d2,h1,'same'); 
g1 = filter(fliplr(h1),1,d2);

M = 12; % length of filter for moving average to get g
h2 = ones(1,M)/M;  % filter for moving average to get g
%g = conv(g1,h2,'same');
g = filter(h2,1,g1);

gth = 0.005; % threshold for  QRS complex detection

suprath = (g> gth); % suprathreshold point
L = 4*(N+M); % window for finding peaks
icnt = 0;
for l =1:length(suprath)
    start = max(l-L,1);
    stop = min(length(suprath),l+L);
    if suprath(l) & (g(l) == max(g(start:stop)))
        icnt = icnt +1;
        peakss(icnt) = l;
    end
end
      

% Plotting Filtered Output
figure(2), clf
subplot(2,1,1)
plot(t,x)
%xlabel('Time [s]')
ylabel('x')
xlim([0 t(end)]);
title('x= ECG')
hold on

subplot(2,1,2)
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

subplot(2,1,1)
hold on
plot(t(suprath),x(suprath),'g.')
plot(t(peakss),x(peakss),'r*')
legend('signal','QRS complex','peaks')

% Extract QRS complexes
iQRS = 0;
for icnt= 2:length(suprath)
    if (suprath(icnt)-suprath(icnt-1))==1
        iQRS = iQRS+1;
        QRS{iQRS} = [];
    end
    if suprath(icnt)
        QRS{iQRS} = [QRS{iQRS} x(icnt)];
    end
end

%% Plot QRS complexes
figure(3), clf
hold on
for iQRS =1 : length(QRS)
    plot((1:length(QRS{iQRS})),QRS{iQRS})
end
xlabel('Time [samples]')
ylabel('Normalised amplitude')
title('Extracted QRS Complexes')

return
%% Analysis: Intuitive measures
maxL = 60;  % max duration in samples
for iQRS =1 : length(QRS)
    iduration(iQRS) = min(length(QRS{iQRS}),maxL); % duration in samples
    duration(iQRS) = min(length(QRS{iQRS}),maxL)*T; % duration in seconds
    height(iQRS) = max(QRS{iQRS}) - min(QRS{iQRS}); 
    midpoint(iQRS) = round(min(length(QRS{iQRS}),maxL)/2);
    centre(iQRS) =  (max(QRS{iQRS}) + min(QRS{iQRS}))/2;
    baseline(iQRS) = (QRS{iQRS}(1)+ QRS{iQRS}(iduration(iQRS)))/2;
    offset(iQRS) = centre(iQRS)-baseline(iQRS);
end

figure(4), clf
subplot(2,2,1)
scatter(offset,height)
xlabel('offset [a.u]')
ylabel('height [a.u.]')
title('ECG features: height vs offset')
subplot(2,2,2)
scatter(offset,duration)
xlabel('offset [a.u]')
ylabel('duration [s]')
title('ECG features: duration vs offset')
return

%% Analysis: PCA
maxL = 34;  % max duration in samples
QRSM = zeros(length(QRS),maxL);
for iQRS =1 : length(QRS)
    iduration = min(length(QRS{iQRS}),maxL);
    QRSM(iQRS,1:iduration) = QRS{iQRS}(1:iduration);
end
                                                         
QRSM = QRSM - mean(QRSM);  % subtract means of each sample time.

C = QRSM'*QRSM;  % correlation/covariance matrix

Q = 2; % number of principal components to calculate
[V D] = eigs(C,Q); % columns of V give principal components = evectors
coef = QRSM*V; % calculate projected  coefficients for PC1 and PC2

figure(4)
subplot(2,2,4)
scatter(coef(:,1),coef(:,2));
xlabel('PCA 1 [a.u.] ')
ylabel('PCA 2 [a.u.] ')
title('PCA projection of QRS complex')

subplot(2,2,3)
plot(V)
xlabel('time [samples]')
ylabel('amplitude [a.u.]')
title('Principal Components')
legend('PC1','PC2')
```
