NewAir = preprocessing('R1',13);
temp = emg('R1',1);
y = temp(:,end);
clear temp
%%
k1 = NewAir(find())
L = 300;
Y = fft(k(1:end-1));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 50*(0:(L/2))/L;
plot(f,P1)
% plot(NewAir(1,:));