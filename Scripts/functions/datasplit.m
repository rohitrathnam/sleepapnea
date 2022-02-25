function [train,test] = datasplit(sub)
%% Setup feature matrix
% feature_matr = [emg_mean emg_rms emg_sd emg_power1 emg_power2 emg_power3 eog_rem eog_sem y];
EMG_feat = [];
EOG_feat = [];
ECG_feat = [];
Resp_feat = [];
EEG_feat = [];
y = [];
disp(sub);
filename = strcat('R',int2str(sub));
temp = emg(filename,1);
y = temp(:,end);
EMG_feat = temp(:,1:end-1);
EOG_feat = eog(filename);
ECG_feat = ecg(filename);
EEG_feat = eeg(filename);
Resp_feat = resp(filename);
% add rest here

%%
feature_matr = [EMG_feat EOG_feat ECG_feat EEG_feat Resp_feat y];
% clear EMG_feat EMG_feat ECG_feat temp Resp_feat EEG_feat;

%% Normalize
k = size(feature_matr);
for i=1:k(2)-1
   feature_matr(:,i) =  (feature_matr(:,i)-min(feature_matr(:,i)))/(max(feature_matr(:,i))-min(feature_matr(:,i)));
end

%% Split test set
ap_ind = find(y==1);
nap_ind = find(y ==0);
m1 = ap_ind(randperm(length(ap_ind),10));
m2 = nap_ind(randperm(length(nap_ind),10));
test = [feature_matr(m1,:); feature_matr(m2,:)];
feature_matr([m1;m2], :) = [];

%% Undersample method
% load('test.mat');
ap_ind = find(feature_matr(:,end)==1);
nap_ind = find(feature_matr(:,end)==0);
[k1,~] = size(nap_ind);
[k2,~] = size(ap_ind);
resamp_nap = datasample(feature_matr(nap_ind,:),k2,1);
set = [feature_matr(ap_ind,:); resamp_nap];
train = set(randperm(size(set, 1)), :);

% %% Oversample method
% ap_ind = find(feature_matr(:,end)==1);
% nap_ind = find(feature_matr(:,end)==0);
% [k1,~] = size(nap_ind);
% [k2,~] = size(ap_ind);
% resamp_ap = datasample(feature_matr(ap_ind,:),k1-k2,1);
% set = [feature_matr; resamp_ap];
% train = set(randperm(size(set, 1)), :);%% Setup feature matrix


% %% ======Checking parameter effectiveness using LASSO Regression============
% [B,FitInfo] = lasso(rset(:,1:61),rset(:,62),'CV',10);
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