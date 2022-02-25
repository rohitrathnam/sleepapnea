function [outputs, labels] = resp(filename)
%% Load and filter
name = strcat(filename,'.mat');
load(name);
ThorRaw = record(9,1:length(stages)*10);
AbdoRaw = record(10,1:length(stages)*10);
AirRaw = record(13,1:length(stages)*10);

% % Equiripple Lowpass filter designed using the FIRPM function.
% 
% % All frequency values are in Hz.
% Fs = 10;  % Sampling Frequency
% 
% Fpass = 0.5;             % Passband Frequency
% Fstop = 1;               % Stopband Frequency
% Dpass = 0.057501127785;  % Passband Ripple
% Dstop = 0.0001;          % Stopband Attenuation
% dens  = 20;              % Density Factor
% 
% % Calculate the order from the parameters using FIRPMORD.
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% 
% % Calculate the coefficients using the FIRPM function.
% b  = firpm(N, Fo, Ao, W, {dens});
% % Hd = dfilt.dffir(b);

%% Low pass filtering
F1 = 0.4; % passband cutoff frequency -> 30 RPM
F2 = 0.6; % stopband cutoff frequency -> 40 RPM
Fs = 10;
d1 = 0.02;
d2 = 0.01;
[n0,f0,a0,w] = firpmord([F1,F2],[1,0],[d1,d2],Fs); % n0 = 174
b = firpm(n0,f0,a0,w);

%%
ThorRes = filtfilt(b,1,ThorRaw);
AbdoRes = filtfilt(b,1,AbdoRaw);
Air = filtfilt(b,1,AirRaw);

% %% Plot
% figure()
% hold on
% plot(ThorRes(1:1000));
% [pks, locs] = findpeaks(ThorRes(1:1000), 'MinPeakDistance', 10);
% scatter(locs, pks);
% hold off

%% Split epochs and annotate
col = 10*30;
row = 1;
tresp_epochs = zeros(length(stages)/30, col+1);
aresp_epochs = zeros(length(stages)/30, col+1);
air_epochs = zeros(length(stages)/30, col+1);
for i = 1:col:length(ThorRes)
    tresp_epochs(row,:) = [ThorRes(1,i:i+col-1) 0];
    aresp_epochs(row,:) = [AbdoRes(1,i:i+col-1) 0];
    air_epochs(row,:) = [Air(1,i:i+col-1) 0];
    row = row+1;
end

Epoch_indices_hypopnea = [];
Epoch_indices_artifact = [];
epoch_idx = 0;
for i = 1:length(Events)
    s = Events(i).EventConcept;
    if strcmp(s,'SignalArtifactEvent') == 1
        idx1 = (Events(i).Start);
        epoch_idx1 = floor(idx1/30)+1;
        idx2 = (Events(i).Start+Events(i).Duration);
        epoch_idx2 = floor(idx2/30)+1;
        if epoch_idx1 ~= epoch_idx2
            if epoch_idx1+1 == epoch_idx2
                epoch_idx = [epoch_idx1 epoch_idx2];
            else
                epoch_idx = linspace(epoch_idx1, epoch_idx2, epoch_idx2-epoch_idx1+1);
            end
        else
            epoch_idx = epoch_idx1;
        end
        Epoch_indices_artifact = [Epoch_indices_artifact epoch_idx];
    elseif strcmp(s,'SDO:HypopneaFinding') == 1
        idx1 = (Events(i).Start);
        epoch_idx1 = floor(idx1/30)+1;
        idx2 = (Events(i).Start+Events(i).Duration);
        epoch_idx2 = floor(idx2/30)+1;
        if epoch_idx1 ~= epoch_idx2
            if epoch_idx1+1 == epoch_idx2
                epoch_idx = [epoch_idx1 epoch_idx2];
            else
                epoch_idx = linspace(epoch_idx1, epoch_idx2, epoch_idx2-epoch_idx1+1);
            end
        else
            epoch_idx = epoch_idx1;
        end
        Epoch_indices_hypopnea = [Epoch_indices_hypopnea epoch_idx];
        %signal_epochs(epoch_idx(1):epoch_idx(length(epoch_idx)), col+1) = 1;
    end
end

for i=1:length(Epoch_indices_hypopnea)
    tresp_epochs(Epoch_indices_hypopnea(i), col+1) = 1;
    aresp_epochs(Epoch_indices_hypopnea(i), col+1) = 1;
    air_epochs(Epoch_indices_hypopnea(i), col+1) = 1;
end
for i=1:length(Epoch_indices_artifact)
    tresp_epochs(Epoch_indices_artifact(i), col+1) = -1;
    aresp_epochs(Epoch_indices_artifact(i), col+1) = -1;
    air_epochs(Epoch_indices_artifact(i), col+1) = -1;
end

%% Find peaks and RR interval
[k,~] = size(tresp_epochs);
mean_RR_int = zeros(k, 3);
num_peaks = zeros(k,3);
for i=1:k
    [pks, locs] = findpeaks(tresp_epochs(i,1:end-1), 'MinPeakDistance', 10);
    num_peaks(i,1) = length(pks);
    for j=2:length(locs)
        mean_RR_int(i,1) = mean_RR_int(i,1) + locs(j) - locs(j-1);
    end
    mean_RR_int(i,1) = mean_RR_int(i,1)/(length(locs)-1);
    
    [pks, locs] = findpeaks(aresp_epochs(i,1:end-1), 'MinPeakDistance', 10);
    num_peaks(i,2) = length(pks);
    if(length(locs)~=1)
        for j=2:length(locs)
            mean_RR_int(i,2) = mean_RR_int(i,2) + locs(j) - locs(j-1);
        end
        mean_RR_int(i,2) = mean_RR_int(i,2)/(length(locs)-1);
    end
    %     hold on
    %     plot(tresp_epochs(i,1:end-1));
    %     scatter(locs, pks);
    %     hold off
    %     close all;
    [pks, locs] = findpeaks(air_epochs(i,1:end-1), 'MinPeakDistance', 10);
    num_peaks(i,3) = length(pks);
    if(length(locs)~=1)
        for j=2:length(locs)
            mean_RR_int(i,3) = mean_RR_int(i,3) + locs(j) - locs(j-1);
        end
        mean_RR_int(i,3) = mean_RR_int(i,3)/(length(locs)-1);
    end
end

% histogram(mean_RR_int_a, [0:1:100]);
% scatter([1:length(mean_RR_int_a)],mean_RR_int_a)

%% Statistical features
std_resp = zeros(k, 3);
skew_resp = zeros(k, 3);
kurt_resp = zeros(k, 3);

for i=1:k
    std_resp(i,1) = std(tresp_epochs(i,1:end-1));
    std_resp(i,2) = std(aresp_epochs(i,1:end-1));
    std_resp(i,3) = std(air_epochs(i,1:end-1));
    skew_resp(i,1) = skewness(tresp_epochs(i,1:end-1));
    skew_resp(i,2) = skewness(aresp_epochs(i,1:end-1));
    skew_resp(i,3) = skewness(air_epochs(i,1:end-1));
    kurt_resp(i,1) = kurtosis(tresp_epochs(i,1:end-1));
    kurt_resp(i,2) = kurtosis(aresp_epochs(i,1:end-1));
    kurt_resp(i,3) = kurtosis(air_epochs(i,1:end-1));
end
outputs = [mean_RR_int num_peaks std_resp skew_resp kurt_resp];
labels = {'mean_resp_intT', 'mean_resp_intA', 'mean_RR_intN', 'num_peaksT', 'num_peaksA', 'num_peaksN', 'std_respT', 'std_respA', 'std_respN', 'skew_respT', 'skew_respA', 'skew_respN', 'kurt_respA', 'kurt_respT', 'kurt_respN'};
end