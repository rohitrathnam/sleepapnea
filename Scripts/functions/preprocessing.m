function signal_epochs = preprocessing(filename,index)
%% Loading data, repeat for each file
% addpath('sleepapnea');
name = strcat(filename,'.mat');
load(name);

%% Segmenting data into 30s epoch

% Extract signals of interest, replace index appropriately
signal = record(index, 1:hdr.samples(index)*length(stages));

% Extracting 30 second epochs
% rows = no. of epochs = length(stages)/30
% columns = no. of samples per epoch + one extra column for apnea state
col = hdr.samples(index)*30;
row = 1;
signal_epochs = zeros(length(stages)/30, col+1);
for i = 1:col:length(signal)
    signal_epochs(row,:) = [signal(1,i:i+col-1) 0];
    row = row+1;
end

%% =====Getting the hypopnea event indices========
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

%%
for i=1:length(Epoch_indices_hypopnea)
    signal_epochs(Epoch_indices_hypopnea(i), col+1) = 1;
end
for i=1:length(Epoch_indices_artifact)
    signal_epochs(Epoch_indices_artifact(i), col+1) = -1;
end

%% ======Getting the hypopnea ECG recordings======
% signal_apnea = [];
% signal_nopnea = []; %Modify where event lasts more than 30s
% for i = 1:size(signal_epochs, 1)
%     if(signal_epochs(i, col+1)) == 1
%         signal_apnea = [signal_apnea; signal_epochs(i,1:col)];
%     else
%         signal_nopnea = [signal_nopnea; signal_epochs(i,1:col)];
%     end
end




