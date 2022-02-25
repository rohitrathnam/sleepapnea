%%
index = 9; % abdominal respiration
signal_epochs = preprocessing('R1',index); 
s = size(signal_epochs);
epochs_nb = s(1);
samples_nb = s(2)-1;
Fs = 10;
%%
figure(1)
hold on;
for i=1:1
    if signal_epochs(i,end) == 0
        plot(signal_epochs(i,:), 'k')
    elseif signal_epochs(i,end) == 1
        plot(signal_epochs(i,:), 'r', 'Linewidth', 1.5)
    end
end
hold off;
%%
fftsignal = fft(signal_epochs(1,:), length(signal_epochs(1,:)));
plot(abs(fftsignal(1:ceil(length(fftsignal)/2))));
%% Low pass filtering
F1 = 0.4; % passband cutoff frequency -> 30 RPM
F2 = 0.6; % stopband cutoff frequency -> 40 RPM
d1 = 0.02;
d2 = 0.01;
[n0,f0,a0,w] = firpmord([F1,F2],[1,0],[d1,d2],Fs); % n0 = 10
f2 = firpm(n0+1,f0,a0,w); % n0+1 to fullfil the specifications
%fvtool(f2,1)

signal_epochs2 = [];
for i=1:epochs_nb
    epoch = filtfilt(f2,1,signal_epochs(i,1:samples_nb-1));
    signal_epochs2(i,:) = [epoch,signal_epochs(i,end)];
end

% plot the original signal and the filtered one
ep = 500;
plot(signal_epochs(ep,:))
hold on;
plot(signal_epochs2(ep,:))
hold off;
%% 
fftnoisy = fft(signal_epochs(1,:), length(signal_epochs(1,:)));
plot(abs(fftnoisy(1:ceil(length(fftnoisy)/2))));
hold on;
fftclean = fft(signal_epochs2(1,:), length(signal_epochs2(1,:)));
plot(abs(fftclean(1:ceil(length(fftclean)/2))), 'r');
hold off;
%% PART II : features selection
%% Maximum value
Max = zeros(1,epochs_nb);
for i=1:epochs_nb
    Max(i) = max(signal_epochs3(i,:));
end

%% plot
for i=1:epochs_nb
    if signal_epochs(i,end) == 0
        plot(1,Max(i), ['b', '.'])
    elseif signal_epochs(i,end) == 1
        plot(1,Max(i), ['r', '_'])
    end
    hold on;
end
hold off;
%% Mean
% square the signal
signal_epochs3 = [];
for i=1:epochs_nb
    signal_epochs3(i,:) = signal_epochs2(i,:).^2;
end

Mn = zeros(1,epochs_nb);
for i=1:epochs_nb
    Mn(i) = mean(signal_epochs3(i,:));
end
%% plot
for i=1:epochs_nb
    if signal_epochs(i,end) == 0
        plot(1,Mn(i), ['b', '.'])
    elseif signal_epochs(i,end) == 1
        plot(1,Mn(i), ['r', '_'])
    end
    hold on;
end
hold off;
%% Median
Md = zeros(1,epochs_nb);
for i=1:epochs_nb
    Md(i) = mean(signal_epochs3(i,1:end-1));
end
plot(signal_epochs3(1,1:end-1));
hold on;
plot(ones(1,samples_nb)*Md(1), 'r');
hold off;
%% plot
for i=1:epochs_nb
    if signal_epochs(i,end) == 0
        plot(1,Md(i), ['b', '.'])
    elseif signal_epochs(i,end) == 1
        plot(1,Md(i), ['r', '_'])
    end
    hold on;
end
hold off;
%% Number of cycles
cycles_nb = zeros(1,samples_nb);
mean_interval_time = zeros(1,samples_nb);
mean_peak_value = zeros(1,samples_nb);
for i=1:epochs_nb
    pks = findpeaks(signal_epochs2(i,:));
    ind = zeros(1,samples_nb);
    for k=1:length(pks)
        ind(k) = find(signal_epochs2(i,:)==pks(k));
    end
    mean_peak_value(i) = mean(pks);
    mean_interval_time(i) = mean(diff(ind));
    cycles_nb(i) = length(pks);
end
%% plot
for i=1:epochs_nb
    if signal_epochs(i,end) == 0
        plot(i,cycles_nb(i), ['b', '.'])
    elseif signal_epochs(i,end) == 1
        plot(i,cycles_nb(i), ['r', '|'])
    end
    hold on;
end
hold off;
%% Power and Energy
power = zeros(1,samples_nb);
energy = zeros(1,samples_nb);
for i=1:epochs_nb
    power(i) = rms(signal_epochs2(i,:))^2;
    energy(i) = norm(signal_epochs2(i,:),2)^2;
end
%% plot
for i=1:epochs_nb
    if signal_epochs(i,end) == 0
        plot(1,energy(i), ['b', '.'])
    elseif signal_epochs(i,end) == 1
        plot(1,energy(i), ['r', '_'])
    end
    hold on;
end
hold off;
%% Features matrix
featMat = [Max; Mn; Md; cycles_nb; mean_interval_time; mean_peak_value; power; energy];
% normalization
for i=1:length(featMat(:,1))
    featMat(i,:) = featMat(i,:)/max(abs(featMat(i,:)));
end
