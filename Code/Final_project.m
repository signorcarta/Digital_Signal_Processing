clear;
clc;
close all;

filename='signal_124.wav';

audioinformazioni=audioinfo(filename); % vectors containing info about the file
nbits=audioinformazioni.BitsPerSample; 
duration=audioinformazioni.Duration; % duration of the audiofile

[y,fs] = audioread(filename); % reads the file

% Sampling period
Tc = 1/fs; % sample period of the audio file (fs=sampling frequency)

% Number of samples on which the fft is computed
N=fs*duration;

% Input signal spectrum____________________________________________________

K = 1/N*fftshift(fft(y(:,1),N)); % Computing the fft and shifting to zero frequency
T = -fs/2:fs/N:fs/2-fs/N;% Creates the frequency axis
figure(1);
plot(T,abs(K));
grid
title('Input signal spectrum');
xlabel('frequency[Hz]');
ylabel('amplitude');

% Carrier frequency detection______________________________________________
data_fft=fft(y);
Y = fftshift(data_fft);
powershift = (abs(Y).^2)/N;
carrier=find(powershift>1000);
fshift = (-N/2:N/2-1)*(fs/N);
f1=fshift(carrier(3));
f2=fshift(carrier(4));

%Survery of the amplitudes at carrier frequency____________________________
amplishift = (abs(Y))./N;
a1=amplishift(carrier(3));
a2=amplishift(carrier(4));

%Extract carrier f1________________________________________________________
delta_f3dB=2;

theta0_f1=2*pi*f1/fs;
delta_theta_3dB_f1=2*pi*delta_f3dB/fs;
delta_f1=delta_theta_3dB_f1/2;
r=1-delta_f1;
c_a1_f1=r*2*cos(theta0_f1);
c_a2_f1=-r^2;
c_b0_f1=(1-r)*2*sin(theta0_f1);
p=r*exp(j*theta0_f1);
z=exp(j*theta0_f1);
%freqz([c_b0_f1],[1 -c_a1_f1 -c_a2_f1])
y_carrierf1 = filter([c_b0_f1],[1 -c_a1_f1 -c_a2_f1], y);

%Extract carrier f2________________________________________________________
theta0_f2=2*pi*f2/fs;
delta_theta_3dB_f2=2*pi*delta_f3dB/fs;
delta_f2=delta_theta_3dB_f2/2;
r=1-delta_f2;
c_a1_f2=r*2*cos(theta0_f2);
c_a2_f2=-r^2;
c_b0_f2=(1-r)*2*sin(theta0_f2);
p=r*exp(j*theta0_f2);
z=exp(j*theta0_f2);
%freqz([c_b0_f2],[1 -c_a1_f2 -c_a2_f2])
y_carrierf2 = filter([c_b0_f2],[1 -c_a1_f2 -c_a2_f2], y);

% Demodulation of carrier f1_______________________________________________
y_demodcarrier1=y.*y_carrierf1;

% Demodulation of carrier f2_______________________________________________
y_demodcarrier2=y.*y_carrierf2;

% Bandpassfilter of carrier1_______________________________________________
[n,fo,ao,w] = firpmord([8000 8500],[1 0],[0.0001 0.001],fs);
b = firpm(n,fo,ao,w);
figure(2);
freqz(b,1);
title('Low-pass FIR filter');

y_demodcarr1_bandpass=filtfilt(b,1,y_demodcarrier1);

% Bandpassfilter of carrier2_______________________________________________

y_demodcarr2_bandpass=filtfilt(b,1,y_demodcarrier2);

% Notch for the origin_____________________________________________________
f0_notch=0;
deltaf3dB_notch=2;
theta0_notch=f0_notch/fs*2*pi;
delta_theta3dB_notch=deltaf3dB_notch/fs*2*pi;
delta_notch=delta_theta3dB_notch/2;
r_notch=1-delta_notch;
b1_notch=-2*cos(theta0_notch);
a1_notch=2*r_notch*cos(theta0_notch);
a2_notch=-r_notch^2;
figure(3);
freqz([1 b1_notch 1],[1 -a1_notch -a2_notch]);
title('High-pass notch IIR filter');

y_final1 = filter([1 b1_notch 1],[1 -a1_notch -a2_notch],y_demodcarr1_bandpass);
y_final2 = filter([1 b1_notch 1],[1 -a1_notch -a2_notch],y_demodcarr2_bandpass);

% Spectrum of the demodulated and filtered signal 1________________________
Y = 1/N*fftshift(fft(y_final1(:,1),N));
F = -fs/2:fs/N:fs/2-fs/N;
figure(4);
plot(F,abs(Y));
grid
xlabel('frequency[Hz]');
ylabel('amplitude');
title('Spectrum of the demodulated and filtered signal 1');

% Spectrum of the demodulated and filtered signal 2________________________

Y = 1/N*fftshift(fft(y_final2(:,1),N));
F = -fs/2:fs/N:fs/2-fs/N;
figure(5);
plot(F,abs(Y));
grid
xlabel('frequency[Hz]');
ylabel('amplitude');
title('Spectrum of the demodulated and filtered signal 2');
% Save ____________________________________________________________________

%player=audioplayer(y_final2,fs,nbits);
%play(player);

% Amplifying the signal louder ____________________________________________
y_final1=y_final1.*10;
y_final2=y_final2.*10;

yfinal=[y_final1(:),y_final2(:)];
audiowrite('signal_output.wav',yfinal,fs);

