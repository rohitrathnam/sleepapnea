function [output,labels] = eog(filename)
EOGR = preprocessing(filename,6);
EOGL = preprocessing(filename,7);
EOG = [EOGR(:,1:end-1)-EOGL(:,1:end-1)];

%% Statistical features
[k,~] = size(EOG);
eog_mean = zeros(k,3);
eog_std = zeros(k,3);
eog_energy = zeros(k,3);
eog_skew = zeros(k,3);
eog_kurt = zeros(k,3);
for sub=1:k
    eog_mean(sub,1) = mean(EOGR(sub,1:end-1));
    eog_mean(sub,2) = mean(EOGL(sub,1:end-1));
    eog_mean(sub,3) = mean(EOG(sub,1:end-1));
    eog_std(sub,1) = std(EOGR(sub,1:end-1));
    eog_std(sub,2) = std(EOGL(sub,1:end-1));
    eog_std(sub,3) = std(EOG(sub,1:end-1));
    eog_skew(sub,1) = skewness(EOGR(sub,1:end-1));
    eog_skew(sub,2) = skewness(EOGL(sub,1:end-1));
    eog_skew(sub,3) = skewness(EOG(sub,1:end-1));
    eog_kurt(sub,1) = kurtosis(EOGR(sub,1:end-1));
    eog_kurt(sub,2) = kurtosis(EOGL(sub,1:end-1));
    eog_kurt(sub,3) = kurtosis(EOG(sub,1:end-1));
    F = fft(EOGR(sub,1:end-1));
    eog_energy(sub,1) = sum(F.*conj(F));
    F = fft(EOGL(sub,1:end-1));
    eog_energy(sub,2) = sum(F.*conj(F));
    F = fft(EOG(sub,1:end-1));
    eog_energy(sub,3) = sum(F.*conj(F));
end

%% Power in bands 0.1-0.3Hz, 0.3-0.5Hz
[k,~] = size(EOG);
eog_01_03 = zeros(k,3);
eog_03_05 = zeros(k,3);
for sub=1:k
    [pxx,f] = pburg(EOGR(sub,1:end-1),10,50);
    eog_01_03(sub, 1) = mean(pxx(2:3));
    eog_03_05(sub, 1) = mean(pxx(4:5));
    [pxx,f] = pburg(EOGL(sub,1:end-1),10,50);
    eog_01_03(sub, 2) = mean(pxx(2:3));
    eog_03_05(sub, 2) = mean(pxx(4:5));
    [pxx,f] = pburg(EOG(sub,:),10,50);
    eog_01_03(sub, 3) = mean(pxx(2:3));
    eog_03_05(sub, 3) = mean(pxx(4:5));
end

% %% Normalize
% eog_rem = (eog_rem - min(eog_rem))/(max(eog_rem) - min(eog_rem));
% eog_sem = (eog_sem - min(eog_sem))/(max(eog_sem) - min(eog_sem));

output = [eog_mean eog_energy eog_std eog_skew eog_kurt eog_01_03 eog_03_05];
labels = {'eog_mean1' 'eog_mean2' 'eog_mean3' 'eog_energyR' 'eog_energyL' 'eog_energyD' 'eog_stdR' 'eog_stdL' 'eog_stdD' 'eog_skewR' 'eog_skewL' 'eog_skewD' 'eog_kurtR' 'eog_kurtL' 'eog_kurtD' 'eog_01_5R' 'eog_01_5L' 'eog_01_5D' 'eog_03_05R' 'eog_03_05L' 'eog_03_05D'};

% %%
% figure()
% subplot(3,1,1)
% plot(EOGR(1,1:end-1));
% subplot(3,1,2)
% plot(EOGL(1,1:end-1));
% subplot(3,1,3)
% plot(EOG(1,1:end-1));

% %% test AR model for best AIC
% fpe = zeros(10, 9);
% for filename=1:10
%     name = strcat(int2str(filename),'.mat');
%     name = strcat('R',name);
%     load(name);
%     signal = record(6, 1:hdr.samples(6)*length(stages));
%     data = iddata(signal',[],0.02);
%     for order=2:10
%         sys = ar(data,order);
%         fpe(filename, order-1) = sys.Report.Fit.FPE;
%         disp([filename, order])
%     end
% end
% 
% 
% %%
% [pxx,f] = pwelch(EOG(5,:), 50);
% data = iddata(EOG(5,:)',[],0.02);
% sys = ar(data,4);
% figure()
% subplot(2,1,1)
% plot(f,pxx)
% subplot(2,1,2)
% spectrum(sys)