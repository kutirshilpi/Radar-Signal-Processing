
% waveform setting

Tau=100e-6;
B=1e6;
Fs=20e6;
t=0:1/Fs:Tau
Fi=B*t/Tau;
xt=exp(j*pi*B*t.*t/Tau)
td=t-Tau;
ht=exp(j*pi*B*td.*td/Tau)

% matched filtering
yt=conv(xt,ht);
X=fft(xt);
X=fftshift(X);
Y=fft(yt);
Y=fftshift(Y);

df = Fs/length(X);
freqvec = -Fs/2+df:df:Fs/2;

%plotting
plot(freqvec/B,abs((X)))
figure
plot(t/Tau,real(xt))
figure
plot(abs((Y)))