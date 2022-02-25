function [outputs, labels] = ox(filename)
sao2_raw = preprocessing(filename,1);
oxstat = preprocessing(filename,14);

% name = strcat(filename,'.mat');
% load(name);
% sao2 = record(1,1:length(stages));
% oxstst = record(14,1:length(stages));
% plot(sao2_raw(1,:))

%% Remove offset and clip; get stat features
[k,~] = size(sao2_raw);
sao2 = zeros(k, 30);
ox_feat = zeros(k, 4);
for i=1:k
    sao2(i,:) = sao2_raw(i,1:end-1)-median(sao2_raw(i,1:end-1));
    for j=1:30
        if sao2(i,j) < -10
            sao2(i,j) = -10;
        end
    end
    ox_feat(i,1) = mean(sao2(i,:));
    ox_feat(i,2) = std(sao2(i,:));
    ox_feat(i,3) = mean(oxstat(i,1:end-1));
    ox_feat(i,4) = std(oxstat(i,1:end-1));
end


%% Return features
outputs = ox_feat;
labels = {'sao2_mean', 'sao2_std', 'ox_mean', 'ox_std'};
end