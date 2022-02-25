function [output,labels] = emg(filename, yflag)
EMG = preprocessing(filename,5);
y = EMG(:,end);
% EOGR = preprocessing('R1',6);
% EOGL = preprocessing('R1',7);
% NewAir = preprocessing('R1',13);

% %% Mean energy bands for apnea and no apnea
% 
% k = size(EMG);
% % average 8 points
% feature = zeros(k(1),5);
% for i=1:k(1)
%     L = 3750;
%     Y = fft(EMG(i,1:end-1));
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = 50*(0:(L/2))/L;
%     %plot(f,P1);
%     %disp(EMG(i,3751));
%     l = 1;
%     for j=1:375:1500
%         feature(j, l) = mean(P1(j:j+375));
%         l = l+1;
%     end
% end

%% Statistical feat
k = size(EMG);
emg_mean = zeros(k(1),1);
% emg_rms = zeros(k(1),1);
emg_sd = zeros(k(1),1);
for i=1:k(1)
    temp = [];
%     for j=1:(k(2)-1)/5
%         temp = [temp rms(EMG(i,(5*(j-1))+1:5*j))];
%     end
    emg_mean(i) = mean(EMG(i,1:end-1));
    emg_sd(i) = std(EMG(i,1:end-1));
end

%% 10s average power (Welch periodogram)

k = size(EMG);
emg_power1 = zeros(k(1),1);
emg_power2 = zeros(k(1),1);
emg_power3 = zeros(k(1),1);

for i=1:k(1)
    emg_power1(i) = mean(pwelch(EMG(i,1:1250)));
    emg_power2(i) = mean(pwelch(EMG(i,1251:2500)));
    emg_power3(i) = mean(pwelch(EMG(i,2501:3750)));
end

% %% Normalisation
% emg_mean = (emg_mean - min(emg_mean))/(max(emg_mean) - min(emg_mean));
% emg_rms = (emg_rms - min(emg_rms))/(max(emg_rms) - min(emg_rms));
% emg_sd = (emg_sd - min(emg_sd))/(max(emg_sd) - min(emg_sd));
% emg_power1 = (emg_power1 - min(emg_power1))/(max(emg_power1) - min(emg_power1));
% emg_power2 = (emg_power2 - min(emg_power2))/(max(emg_power2) - min(emg_power2));
% emg_power3 = (emg_power3 - min(emg_power3))/(max(emg_power3) - min(emg_power3));

% %% Set output
% if(yflag)
%     output = [emg_mean emg_rms emg_sd emg_power1 emg_power2 emg_power3 y];
% else
%     output = [emg_mean emg_rms emg_sd emg_power1 emg_power2 emg_power3];
% end

%% Set output
if(yflag)
    output = [emg_mean emg_sd emg_power1 emg_power2 emg_power3 y];
else
    output = [emg_mean emg_sd emg_power1 emg_power2 emg_power3];
end

labels = {'emg_mean', 'emg_sd', 'emg_power1', 'emg_power2', 'emg_power3'};
end