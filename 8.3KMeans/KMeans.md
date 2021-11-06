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

% down sample PSD by a factor of 10
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

%% Calculate principle components
%SdB = SdB - mean(SdB(:));                                                  % INCORRECT: remove mean for calculating variance about mean
Sav = mean(SdB,2);                                                          % Calculate mean value of each dimension
SdB = SdB - repmat(Sav,1,size(SdB,2)); 
C = SdB*SdB'/size(SdB,2);  % correlation/covariance matrix

Q = 10; % number of principle components to calculate
[V D] = eigs(C,Q); % columns of V give principle components = evectors
                   %
% change sign of evectors if the are larger at high than low frequecny
for iv = 1:Q
    V(:,iv) = sign(V(1,iv)-V(end,iv))*V(:,iv);
end

%% Plot principle component 
figure(1)
subplot(2,1,2)
q = 3;   % number of component to use in plot 
plot(f, V(:,1:q),'LineWidth', 2)
xlabel('Freq [Hz]')
ylabel('PSD [dB/Hz]')
title('Principle components of vowels')
legend('1st','2nd','3rd')

%% Project  vowel 
coef = SdB'*V;
q = 2;   % number of components to use in clustering

figure(2)
cmap = jet(10);
if q ==1
    scatter(coef(:,1),0.1*randn(size(coef,1),1),[],cmap(idn,:));
    ylim([-1 1])
    xlabel('Principle Component 1 [dB/Hz] ')
elseif q==2
    scatter(coef(:,1),coef(:,2),[],cmap(idn,:));
    xlabel('Principle Component 1 [dB/Hz] ')
    ylabel('Principle Component 2 [dB/Hz] ')
elseif q==3
    scatter3(coef(:,1),coef(:,2),coef(:,3),[],cmap(idn,:));
    xlabel('Principle Component 1 [dB/Hz] ')
    ylabel('Principle Component 2 [dB/Hz] ')
    zlabel('Principle Component 3 [dB/Hz] ')
else
    scatter(coef(:,1),coef(:,2),[],cmap(idn,:));
    xlabel('Principle Component 1 [dB/Hz] ')
    ylabel('Principle Component 2 [dB/Hz] ')
end
title('True PC projections of vowels')

    
 %% Clustering 
 K = 9;  % number of clusters
 [kdx, cnt, sumd] = kmeans(coef(:,1:q), K,'Start','plus');
 
figure(3)
cmap = jet(K);
if q ==1
    scatter(coef(:,1),0.1*randn(size(coef,1),1),[],cmap(kdx,:));
       ylim([-1 1])
    xlabel('Principle Component 1 [dB/Hz] ')
elseif q==2
    scatter(coef(:,1),coef(:,2),[],cmap(kdx,:));
    xlabel('Principle Component 1 [dB/Hz] ')
    ylabel('Principle Component 2 [dB/Hz] ')
elseif q==3
    scatter3(coef(:,1),coef(:,2),coef(:,3),[],cmap(kdx,:));
    xlabel('Principle Component 1 [dB/Hz] ')
    ylabel('Principle Component 2 [dB/Hz] ')
    zlabel('Principle Component 3 [dB/Hz] ')
else
    scatter(coef(:,1),coef(:,2),[],cmap(kdx,:));
    xlabel('Principle Component 1 [dB/Hz] ')
    ylabel('Principle Component 2 [dB/Hz] ')
end
title('K-means clustering of PC projections of vowels')

%return
% %% Investigate Number of Clusters K
% figure(4), clf
% for K =2:16
%     cmap = jet(K);
%     
%     [kdx, cnt, sumd] = kmeans(coef(:,1:q), K,'Start','plus'); % sumd give with cluster sum of squares
%     
%     subplot(4,4,K)
%     hold on
%     scatter(coef(:,1),coef(:,2),[],cmap(kdx,:));
%     %plot(cnt(:,1),cnt(:,2),'g*')
%     xlabel('PC1 [dB/Hz] ')
%     ylabel('PC2 [dB/Hz] ')
%     title(['K=' num2str(K)])
% end
```
