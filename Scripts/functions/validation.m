% Add full folder to path
load('workspace.mat');

%% Validation for normal models, trained with the rset variable
y = val_set(:,end);
ap_ind = find(y==1);
nap_ind = find(y==0);

[k1,~] = size(nap_ind);
[k2,~] = size(ap_ind);
resamp_nap = datasample(val_set(nap_ind,:),k2,1);
set = [val_set(ap_ind,:); resamp_nap];

k = size(val_set);
y = svm.predictFcn(set(:,1:k(2)-1));
y_true = set(:,k(2));

C = confusionmat(y, y_true)
confusionchart(C);

%%
y = val_set(:,end);
ap_ind = find(y==1);
nap_ind = find(y==0);

[k1,~] = size(nap_ind);
[k2,~] = size(ap_ind);
resamp_nap = datasample(val_set(nap_ind,:),k2,1);
set = [val_set(ap_ind,:); resamp_nap];

k = size(val_set);
y = bagged.predictFcn(set(:,1:k(2)-1));
y_true = set(:,k(2));

C = confusionmat(y, y_true)
confusionchart(C);

%% Sparse version validation, trained with sparse_set
y = val_set(:,end);
ap_ind = find(y==1);
nap_ind = find(y==0);

[k1,~] = size(nap_ind);
[k2,~] = size(ap_ind);
resamp_nap = datasample(val_set(nap_ind,:),k2,1);
set = [val_set(ap_ind,:); resamp_nap];

k = size(set);
y = svmsparse.predictFcn(set(:, B(:,idxLambda1SE)~=0));
y_true = set(:,k(2));

C = confusionmat(y, y_true)
confusionchart(C);

%%
y = val_set(:,end);
ap_ind = find(y==1);
nap_ind = find(y==0);

[k1,~] = size(nap_ind);
[k2,~] = size(ap_ind);
resamp_nap = datasample(val_set(nap_ind,:),k2,1);
set = [val_set(ap_ind,:); resamp_nap];

k = size(set);
y = baggedsparse.predictFcn(set(:, B(:,idxLambda1SE)~=0));
y_true = set(:,k(2));

C = confusionmat(y, y_true)
confusionchart(C);