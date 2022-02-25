%% Setup feature matrix
val_split = 0.15; % validation split ratio

EMG_feat = [];
EOG_feat = [];
ECG_feat = [];
Resp_feat = [];
EEG_feat = [];
Ox_feat = [];
train_set_all = [];
val_set = [];
y = [];
% for i=1:10
%     disp(i);
%     filename = strcat('R',int2str(i));
%     temp = emg(filename,1);
%     y = [y; temp(:,end)];
%     EMG_feat = [EMG_feat; temp(:,1:end-1)];
%     EOG_feat = [EOG_feat; eog(filename)];
%     ECG_feat = [ECG_feat; ecg(filename)];
%     EEG_feat = [EEG_feat; eeg(filename)];
%     Resp_feat = [Resp_feat; resp(filename)];
%     Ox_feat = [Ox_feat; ox(filename)];
%     % add rest here
% end

for i=1:10
    disp(i);
    filename = strcat('R',int2str(i));
    [temp, l1] = emg(filename,1);
    y = temp(:,end);
    EMG_feat = temp(:,1:end-1);
    [EOG_feat, l2] = eog(filename);
    [ECG_feat, l3] = ecg(filename);
    [EEG_feat, l4] = eeg(filename);
    [Resp_feat, l5] = resp(filename);
    [Ox_feat, l6] = ox(filename);
    train_set = [EMG_feat EOG_feat ECG_feat EEG_feat Resp_feat Ox_feat y];
    train_set(find(y==-1), :) = []; %remove artefact epochs
    
    % Normalize
    k = size(train_set);
    for i=1:k(2)-1
        m1 = min(train_set(:,i));
        m2 = max(train_set(:,i));
        for ii=1:k(1)
            train_set(ii,i) =  (train_set(ii,i)-m1)/(m2-m1);
        end
    end
    
    ap_ind = find(y==1);
    nap_ind = find(y==0);
    
    [k1,~] = size(nap_ind);
    [k2,~] = size(ap_ind);
    nap_count = round(k1*val_split);
    ap_count = round(k2*val_split);
    [vals, idx] = datasample(train_set, nap_count,1);
    train_set(idx, :) = [];
    val_set = [val_set; vals];
    [vals, idx] = datasample(train_set, ap_count,1);
    train_set(idx, :) = [];
    val_set = [val_set; vals];

    train_set_all = [train_set_all; train_set];
end

labels = [l1 l2 l3 l4 l5 l6];

%% Random undersample method
% load('test.mat');
y = train_set_all(:,end);
ap_ind = find(y==1);
nap_ind = find(y==0);

[k1,~] = size(nap_ind);
[k2,~] = size(ap_ind);
resamp_nap = datasample(train_set_all(nap_ind,:),k2,1);
set = [train_set_all(ap_ind,:); resamp_nap];
rset = set(randperm(size(set, 1)), :);

% %% Oversample method
% ap_ind = find(y==1);
% nap_ind = find(y ==0);
% art_ind = find(y==-1);
% [k1,~] = size(nap_ind);
% [k2,~] = size(ap_ind);
% resamp_ap = datasample(train_set_all(ap_ind,:),k1-k2,1);
% set = [train_set_all(ap_ind,:); train_set_all(nap_ind,:); resamp_ap];
% rset = set(randperm(size(set, 1)), :);
% 
% %%
% 
% rset = train_set_all(randperm(size(train_set_all,1)), :);


%% ======Checking parameter effectiveness using LASSO Regression============

[B,FitInfo] = lasso(train_set_all(:,1:end -1),train_set_all(:,end),'CV',10,'PredictorNames',labels);
%lassoPlot(B,FitInfo,'PlotType','CV');

% Display the variables in the model that corresponds to the minimum 
% cross-validated mean squared error (MSE).
idxLambdaMinMSE = FitInfo.IndexMinMSE;
minMSEModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinMSE)==0);
% Gives MHR, RMSSD and SDNN as the output. If you just drop the MRR, the
% accuracy of the model is not affected. It reamins the same

% Display the variables in the sparsest model within one standard error 
% of the minimum MSE.
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0);
% Gives MRR and MHR as the output. If we just use these 2 features, the 
% accuracy did improve.

%%
train_set_sparse = train_set_all(:, B(:,idxLambda1SE)~=0);
train_set = [train_set_sparse train_set_all(:,end)]; 

y = train_set(:,end);
ap_ind = find(y==1);
nap_ind = find(y==0);

[k1,~] = size(nap_ind);
[k2,~] = size(ap_ind);
resamp_nap = datasample(train_set(nap_ind,:),k2,1);
set = [train_set(ap_ind,:); resamp_nap];
sparse_set = set(randperm(size(set, 1)), :);
