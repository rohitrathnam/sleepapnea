%%
index = 3; % EEG
signal_epochs = preprocessing(filename,index); 
s = size(signal_epochs);
epochs_nb = s(1);
samples_nb = s(2)-1;
Fs = 125;
%%
figure(1)
hold on;
for i=1:epochs_nb
    if signal_epochs(i,end) == 0
        plot(signal_epochs(i,:), 'k')
    elseif signal_epochs(i,end) == 1
        plot(signal_epochs(i,:), 'r', 'Linewidth', 1.5)
    end
end
hold off;
%% LOW PASS FILTER
F1 = 30; % passband cutoff frequency
F2 = 40; % stopband cutoff frequency
d1 = 0.02;
d2 = 0.01;
[n0,f0,a0,w] = firpmord([F1,F2],[1,0],[d1,d2],Fs); % n0 = 10
f2 = firpm(n0+1,f0,a0,w); % n0+1 to fullfil the specifications
%delay = mean(grpdelay(f2));
%fvtool(f2,1)
signal_epochs2 = zeros(epochs_nb, samples_nb);
for i=1:epochs_nb
    epoch = filtfilt(f2,1,signal_epochs(i,1:samples_nb-1));
    signal_epochs2(i,:) = [epoch,signal_epochs(i,end)];
end

% plot the original signal and the filtered one
plot(signal_epochs(25,:))
hold on;
plot(signal_epochs2(25,:))
hold off;
%% 
fftnoisy = fft(signal_epochs(1,:), length(signal_epochs(1,:)));
plot(abs(fftnoisy(1:ceil(length(fftnoisy)/2))));
hold on;
fftclean = fft(signal_epochs2(1,:), length(signal_epochs2(1,:)));
plot(abs(fftclean(1:ceil(length(fftclean)/2))), 'r');
hold off;
%%
plot(periodogram(signal_epochs2(1,:)))
% phi = K/N fs
% K = N*phi/fs
%% Separate the signal in different bands (periodogram)
delta = [];
theta = [];
alpha = [];
beta = [];
gamma = []; 
for i=1:epochs_nb
    [pxx,f] = pwelch(signal_epochs2(i,:),[],[],[],Fs);
    f1 = find(f>0&f<=4);
    delta(:,i) = mean(pxx(f1(1):f1(end)));
    f2 = find(f>4&f<=8);
    theta(:,i) = mean(pxx(f2(1):f2(end)));
    f3 = find(f>8&f<=14);
    alpha(:,i) = mean(pxx(f3(1):f3(end)));
    f4 = find(f>14&f<=30);
    beta(:,i) = mean(pxx(f4(1):f4(end)));
    f5 = find(f>30);
    gamma(:,i) = mean(pxx(f5(1):f5(end)));
end
%% PART II : features selection
%% Mean
Mn = zeros(1,epochs_nb);
for i=1:epochs_nb
    Mn(i) = mean(signal_epochs2(i,:));
end
plot(signal_epochs2(1,:));
hold on;
plot(ones(1,samples_nb)*Mn(1), 'r');
hold off;
%% Median
Md = zeros(1,epochs_nb);
for i=1:epochs_nb
    Md(i) = mean(signal_epochs2(i,1:end-1));
end
plot(signal_epochs2(1,1:end-1));
hold on;
plot(ones(1,samples_nb)*Md(1), 'r');
hold off;
%% Band energy ratio
Rdt = delta./theta;
Rda = delta./alpha;
Rdb = delta./beta;
Rdg = delta./gamma;
Rta = theta./alpha;
Rtb = theta./beta;
Rtg = theta./gamma;
Rab = alpha./beta;
Rag = alpha./gamma;
Rbg = beta./gamma;
%% Features matrix
%featMat = [Mn; Md; Edelta; Etheta; Ealpha; Ebeta; Egamma; Rdt; Rda; Rdb; Rdg; Rta; Rtb; Rtg; Rab; Rag; Rbg; Pdelta; Ptheta; Palpha; Pbeta; Pgamma; Vdelta; Vtheta; Valpha; Vbeta; Vgamma];
featMat = [Mn; Md; Rdt; Rda; Rdb; Rdg; Rta; Rtb; Rtg; Rab; Rag; Rbg; delta; theta; alpha; beta; gamma];

% Normalization of the features
for i=1:length(featMat(:,1))
    featMat(i,:) = featMat(i,:)/max(abs(featMat(i,:)));
end

