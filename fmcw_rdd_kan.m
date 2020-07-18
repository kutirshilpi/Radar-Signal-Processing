% Written By : Ashikur Rahman
% FMCW Range Doppler Processing
clear
close all
BW = 0.2e9; % Sweep Bandwidth
T = 35.6E-6; % Pulse repeatition period
fc = 79E9; % Center frequency


c = 3E8; % Speed of light
n = 512; % Number of chirp
NFFT_r = 1024; % range FFT size
M = 1024; % Number of samples for each chirp in data acquisition
Tx = 3; % Transmitter channels
Rx = 4; % Number of recieve channels
al = BW/T; % Slope of ramp
R = 60; % Test range
Fs = 20E6;
ts = 0: 1/Fs: T;% Time vector of range sampling
fnt = 12; % font size for plots
v = 10; % Velocity 
v2=5;
xm = 10;
NFFT_d = 256;
K = 256;
% phase_noise_freq = [ 1e3, 10e3, 100e3:100e3:900e3, 1e6:1e6:, 9.99e6 ]; % Offset From Carrier
% phase_noise_power = [ -84, -100, -96, -109, -122 ]; % Phase Noise power
load phaseNoiseMeas.mat;
phase_noise_freq = phaseNoiseMeas(1:1700,1); 
% phase_noise_power = 0*phaseNoiseMeas(1:1700,2)-200;
phase_noise_power = phaseNoiseMeas(1:1700,2);


% theory
% delta_Phi = 10*log10(2 * (sin(pi*phase_noise_freq*2*R/c)).^2 .* Fs/NFFT_r) + phase_noise_power;
delta_Phi = 10*log10(2 * (sin(pi*phase_noise_freq*2*R/c)).^2) + phase_noise_power;
% delta_Phi = delta_Phi - 10*log10(sum(rangeWin)/M1);
r=.5
figure;
plot(phase_noise_freq,[delta_Phi phase_noise_power])
for n=1:K

    fc = 79E9;
%       fc=fc+randi([1,K])*1e6;
fc=fc+n*4e6;

x1(:,n)=exp(1i*(2*pi*(2*al*R*ts/c+2*fc*v*n*T/c)+4*pi*fc*R/c)); % 2D chirp vector time domain

%      x1(:,n)=exp(1i*(2*pi*(2*al*R*ts/c+2*fc*v*ts/c)+4*pi*fc*R/c)); % 2D chirp vector time domain
    
%      x2(:,n)=exp(1i*(2*pi*(2*al*(R+r)*ts/c+2*fc*v*ts/c)+4*pi*fc*(R+r)/c));
    
      x2(:,n)=exp(1i*(2*pi*(2*al*(R+r)*ts/c+2*fc*v2*n*T/c)+4*pi*fc*(R+r)/c));
     
     x1(:,n)= x1(:,n)+x2(:,n);
     Sout = add_phase_noise( x1(:,n), Fs, phase_noise_freq, phase_noise_power,0 );
%     Sout = add_phase_noise( x1(:,n), Fs, phase_noise_freq, delta_Phi,0 );
   
%     n_shift = round( 2*(R+v*ts)/c * Fs );
%     n_shift = round( 2*R/c * 20e6 );
%     Sout_down = Sout(n_shift+1:end) .* conj(Sout(1:end-n_shift));
%     x(:,n) = x1(n_shift+1:end,n) .* Sout_down;
%     
    clear x;
    x(:,n)=Sout;
%     x(:,n)=x1(:,n);
    rangeWin = chebwin( size(x,1), 100 );
    y(:,n)= fft( x(:,n).*rangeWin, NFFT_r ); % FFT of each chirp, range FFT
end

%%doppler FFT first
yd=zeros(size(y));
M1 = size(y,1);
dopWin = chebwin( NFFT_d, 100 );
for md = 1 : M1
%   x(md,:)=5/2*cos(2*pi*(2*al*R*ts/c+2*fc*v*n*T/c)+4*pi*fc*R/c)+xm;
    yd(md,:) = fft(y(md,:).*dopWin.',NFFT_d); % Doppler FFT
end

% theory
delta_Phi = 10*log10(2 * (sin(pi*phase_noise_freq*2*R/c)).^2 .* Fs/NFFT_r) + phase_noise_power;
delta_Phi = delta_Phi - 10*log10(sum(rangeWin)/M1);

% Plot
figure;
plot(phase_noise_freq, delta_Phi);
xlabel('Offset Frequency (Hz)');
ylabel('dBc');
ylim([-100, 0]);

figure; 
surf(db(yd)); 
shading flat;

figure; 
colormap jet;
    imagesc(db(yd'));
%     caxis([caxis_lower_limit caxis_upper_limit]);




figure;
hold on;
plot( db(max(abs(yd),[],2 )) );
plot( db( mean( abs(y), 2 ) ) );
xlabel('range bin index');
ylabel('magnitude(dB)');
legend('MPRB','Averaged range FFT');

figure;
plot( db(yd(230,:)) );
xlabel('Doppler bin index');
ylabel('magnitude(dB)');
title('Doppler cut');


figure;
plot(phase_noise_freq,[delta_Phi phase_noise_power])

% ///////////////////////////////////////////
% 
% figure
% surf(abs(x1))
% set(gca,'FontSize',fnt)
% 
% figure
% surf(abs(y)) % plot range FFT
% xlabel('Chirps')
% ylabel('Range FFT samples')
% zlabel('Amplitude')
% title('Range FFT stage')
% set(gca,'FontSize',fnt)
