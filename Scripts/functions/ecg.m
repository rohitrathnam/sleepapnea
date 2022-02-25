function [output, labels] = ecg(filename)
% Get data
ECG_epochs = preprocessing(filename, 4);

%% Checking the heart rate variability for the apnea and non apnea episodes
load('num.mat');%Load filter coefficients
% Filter Specifications
% FIR Bandpass filter
% Fs = 125Hz
% Fstop1 = 8Hz
% Fpass1 = 10Hz
% Fpass2 = 25Hz
% Fstop2 = 27Hz
% Stopband attenuation = 60dB
% Passband attenuation = 1dB
% Filter order = 494
% Used these specifications because a lot of papers had used these
% parameters for a bandpass filter. Some had 10Hz to 45Hz while others had 10 to 25HZ
%%
for i = 1:size(ECG_epochs,1)
    ECG_epochs(i,1:3750) = ECG_epochs(i,1:3750) - mean(ECG_epochs(i,1:3750));%Removing the DC components
    ECG_epochs(i,1:3750) = filter(num,1,ECG_epochs(i,1:3750));%Bandpass filtering
end

ECG_apnea = [];
ECG_nopnea = [];
for i = 1:size(ECG_epochs,1)
    if ECG_epochs(i,3751) == 1
        ECG_apnea = [ECG_apnea; ECG_epochs(i,1:end-1)];
    else
        ECG_nopnea = [ECG_nopnea; ECG_epochs(i,1:end-1)];
    end
end
% A delay is introduced by filtering the signal.
% We need to compensate for that
delay1 = mean(grpdelay(num,1,500)); %where num = filter coefficient
% delay1 = 247 samples
% 1 = denominator and 500 = Sampling frequency
sf = ECG_epochs;
sf(:,1:delay1) = []; % Dropping the first 247 samples
h = 1/125; %where h is th step size and fs = 125Hz
Der = diff(sf(:,1:end - 1),1,2)/h; %Differentiating along the rows and dropping
%the last column since it contains 1 and 0s for the labels
Der_sqrd = Der.^2;
% Moving Integrator =======================================================
window_len = 18; %Setting window length for convolution
t = max(Der_sqrd, [], 2); %Return the largest element across each row

delay2 = floor(window_len/2); % delay in convolution is half the delay length
N = length(Der_sqrd);

%% Skip
% Integrated_Signal = [];
% for i = 1:size(Der_sqrd,1)
%     rect = ones(1 ,window_len)./t(i);
%     Int_sig = conv(Der_sqrd(i,:) ,rect);
%     Int_sig = Int_sig(delay2: N);
%     Int_sig = Int_sig./max(Int_sig);% Normalizing the signal
%     Integrated_Signal = [Int_sig; Integrated_Signal];
% end

%% Finding peaks===========================================================
% Sig_index = [];
% for i = 1:size(Integrated_Signal,1)
%     m1 = max(Integrated_Signal(i,:));
%     m2 = mean(Integrated_Signal(i,:));
%     New_sig = (Integrated_Signal(i,:)>(m1*m2)); %Setting threshold as mean x max value of signal
%     leading = find(diff(New_sig)==1); %finding the indices of the leading edge
%     trailing = find(diff(New_sig)==-1);%finding the indices of the trailing edge
%     if leading(1)<trailing(1)
%         leading = leading(1:length(trailing));% As we want a pair of leading and trailing edges
%         %We don't consider the last signals where the signal has detected the leading edge but not the trailing edge
%         % Extrating the signal indices with the maximum values corresponding to the R Waves
%     else
%          trailing = trailing(2:length(trailing));
%          leading = leading(1:length(trailing));
%     end
%     SI = zeros(1,300);
%     for x=1:length(leading)
%         [R_value, R_loc] = max(sf(i, leading(x):trailing(x)));
%         SI(x) = R_loc+leading(x)-1;
%     end
%     Sig_index = [Sig_index; SI];
% end
% 
% hold on
% plot(sf(516,:));
% scatter(Sig_index(516,:), 0.06*ones(300,1));
% hold off;

% Using findpeaks function=================================================
% Extracting time domain features
Sig_index = [];
Peak_heights = [];
RR_dist = [];

MRR = []; % Mean of RR intervals
MHR = []; % Mean heart rate 
RMSSD = []; % Root mean square of differences between adjacent RR intervals
SDNN = []; % Standard deviation of rr intervals
NN50 = []; % Number of adjacent RR intervals exceeding 50 milliseconds
PNN50 = []; % NN50/ Number of RR intervals
len_pow = floor(size(sf,2)/3);
ECG_POW1 = [];
ECG_POW2 = [];
ECG_POW3 = [];

for i = 1:size(ECG_epochs,1)
    SI = zeros(1,60);
    PH = zeros(1,60);
    RR = zeros(1,60);
    [pks, loc] = findpeaks(sf(i,:),'MinPeakDistance',60, 'MinPeakHeight', (0.04 + mean(sf(i,:))));
    p_len = length(pks);
    Dist = diff(loc);
    
    MRR = [MRR; mean(Dist)];
    MHR = [MHR; 2*p_len];
    RMSSD = [RMSSD; sqrt(sum(Dist.*Dist)/length(Dist))];
    SDNN = [SDNN; std(Dist)];
    NN50 = [NN50; sum(Dist<=13)];  
    PNN50 = [PNN50; sum(Dist<=13)/length(Dist)]; % NN50/ Number of RR intervals
    ECG_POW1 = [ECG_POW1; mean(pwelch(sf(i,1:len_pow)))];
    ECG_POW2 = [ECG_POW2; mean(pwelch(sf(i,len_pow+1:2*len_pow)))];
    ECG_POW3 = [ECG_POW3; mean(pwelch(sf(i,2*len_pow+1:end)))];
    
    for x=1:length(loc)
        SI(x) = loc(x);
    end
    
    for x=1:p_len
        PH(x) = pks(x);
    end
    
    for x=1:length(Dist)
        RR(x) = Dist(x);
    end
    Sig_index = [Sig_index; SI];
    Peak_heights = [Peak_heights; PH];
    RR_dist = [RR_dist; RR];
end

% % Min max normalization
% for i = 1:length(MRR);
%     MRR(i) = (MRR(i) - min(MRR))/(max(MRR) - min(MRR));
%     MHR(i) = (MHR(i) - min(MHR))/(max(MHR) - min(MHR));
%     RMSSD(i) = (RMSSD(i) - min(RMSSD))/(max(RMSSD) - min(RMSSD));
%     SDNN(i) = (SDNN(i) - min(SDNN))/(max(SDNN) - min(SDNN));
%     NN50(i) = (NN50(i) - min(NN50))/(max(NN50) - min(NN50));
%     PNN50(i) = (PNN50(i) - min(PNN50))/(max(PNN50) - min(PNN50));
% end

output = [MRR MHR RMSSD SDNN ECG_POW1 ECG_POW2 ECG_POW3];
labels = {'MRR', 'MHR', 'RMSSD', 'SDNN', 'ECG_POW1', 'ECG_POW2', 'ECG_POW3'};

%% Fix imbalance
% Feat_mat_1 = [];
% Feat_mat_2 = [];
% thres = sum(ECG_epochs(:,end)); % threshold = 84 in this case. 
% count = 0;
% 
% %For patient 1 we have 1000 non apnea and 84 apnea features
% % Extracting the  84 Apnea features 
% for i = 1:size(ECG_epochs, 1)
%     if(Feat_mat(i,end)==1) && count <= thres
%         Feat_mat_1 = [Feat_mat_1; Feat_mat(i, :)];
%          count = count + 1;
%     end
% end
% 
% % Extracting 84 Non-Apnea features 
% count = 0;
% for i = 1:size(ECG_epochs, 1)
%     if(Feat_mat(i,end)==0) && count <= thres
%         Feat_mat_2 = [Feat_mat_2; Feat_mat(i, :)];
%         count = count + 1;
%     end
% end
% 
% new_mat = [Feat_mat_1; Feat_mat_2];
% random_x = new_mat(randperm(size(new_mat, 1)), :);
% 
% hold on
% plot(sf(345,:));
% scatter(Sig_index(345,:), 0.06*ones(60,1));
% hold off;
% 
% 
% %% ======Checking parameter effectiveness using LASSO Regression============
% [B,FitInfo] = lasso(random_x(:,1:4),random_x(:,5),'CV',10,'PredictorNames',{'MRR','MHR','RMSSD','SDNN'});
% %lassoPlot(B,FitInfo,'PlotType','CV');
% 
% % Display the variables in the model that corresponds to the minimum 
% % cross-validated mean squared error (MSE).
% idxLambdaMinMSE = FitInfo.IndexMinMSE;
% minMSEModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinMSE)~=0)
% % Gives MHR, RMSSD and SDNN as the output. If you just drop the MRR, the
% % accuracy of the model is not affected. It reamins the same
% 
% % Display the variables in the sparsest model within one standard error 
% % of the minimum MSE.
% idxLambda1SE = FitInfo.Index1SE;
% sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
% % Gives MRR and MHR as the output. If we just use these 2 features, the 
% % accuracy did improve.

end
