```
clear all
%% Load vowel PSD envelops 
% env contains PSD envelopes for:
%   10 different vowels
%   32 repeats of each vowel (each with a different "speaker")
% as 10*32 rows of envelope data.
% str.vowels contains name of the vowels
% idn contains the identity 1-10 of the vowel in rows of env
% f_Hz is the frequency sampling for the spectral envelope

load 'VowelData.mat'

% downsample PSD by a factor of 10
SdB = env(:,1:10:end);
SdB = SdB';
f = f_Hz(1:10:end);

% plot all the spectral envelopes 
figure(1), clf
subplot(2,1,1)
plot(f, SdB)
xlabel('Freq [Hz]')
ylabel('PSD [dB/Hz]')
title('PSD of 10 vowels, 32 repeats')

%return
%% Calculate principal components
%SdB = (SdB' - mean(SdB'))';  % corrected: now removes mean across repeats for each time sample calculating variance about
C = SdB*SdB'/size(SdB,2);  % correlation/covariance matrix

Q = 10; % number of principal components to calculate
[V D] = eigs(C,Q); % columns of V give principal components = evectors
                   %
% change sign of evectors if the are larger at high than low frequecny
for iv = 1:Q
    V(:,iv) = sign(V(1,iv)-V(end,iv))*V(:,iv);
end

%% Plot principal component 
figure(1)
subplot(2,1,2)
q = 3;   % number of component to use in plot 
plot(f, V(:,1:q),'LineWidth', 2)
xlabel('Freq [Hz]')
ylabel('PSD [dB/Hz]')
title('Principle components of vowels')
legend('1st','2nd','3rd')

%% Project all the vowel onto the principal components & plot in PC space
coef = SdB'*V;

figure(2)
cmap = jet(10);
scatter(coef(:,1),coef(:,2),[],cmap(idn,:)); % 2D colour coded
%scatter(coef(:,1),coef(:,2));  % 2D non-colour coded
%scatter3(coef(:,1),coef(:,2),coef(:,3),[],cmap(idn,:)); % 3D colour coded
%scatter(coef(:,1),randn(size(coef,1),1),[],cmap(idn,:)); %1D colour coded
xlabel('Principal Component 1 [dB/Hz] ')
ylabel('Principal Component 2 [dB/Hz] ')
title('Principal component projections of vowels')

%% Reconstructions of PSD envelopes from projections
figure(3), clf
q = 3;   % number of component to use in reconstruction 
v = V(:,1:q);
cf = coef(:,1:q);
recon = v*cf';
for iv =1:9
    subplot(3,3,iv)
    hold on
    col = 32*iv;
    plot(f,  SdB(:, col),'k','LineWidth', 2)
    plot(f,  recon(:, col),'r:','LineWidth', 2)
    title(str.vowels(iv,:))
end
subplot(3,3,1)
xlabel('Freq [Hz]')
ylabel('PSD [dB/Hz]')
legend('original','recontructed')
```
    


