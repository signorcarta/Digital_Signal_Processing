clear;
clc;
close all;

file='signal_124.wav';

vect_info=audioinfo(file); % vectors containing info about the file 
duration=vect_info.Duration; % duration of the audiofile
[y,fs] = audioread(file); % reads the file
Tc = 1/fs; % sampling period 
N=fs*duration; % # samples on which the fft is computed

%% Input signal spectrum___________________________________________________
K = 1/N*fftshift(fft(y(:,1),N)); % Computing fft and shifting to 0-frequency
T = -fs/2:fs/N:fs/2-fs/N;% Creating the frequency axis

% Plotting input signal spectrum
figure(1);
plot(T,abs(K));
grid
title('Input signal spectrum');
xlabel('frequency[Hz]');
ylabel('amplitude');

%% Carrier frequency detection_____________________________________________
data_fft=fft(y); 
Y = fftshift(data_fft);
powershift = (abs(Y).^2)/N;
carrier=find(powershift>1000);
f_shift = (-N/2:N/2-1)*(fs/N);
f1 = f_shift(carrier(3));
f2 = f_shift(carrier(4));

%% Carrier frequency amplitudes____________________________________________
ampli_shift = (abs(Y))./N;
A1=ampli_shift(carrier(3)); 
A2=ampli_shift(carrier(4)); 

%% Extract carriers________________________________________________________
% Extract carrier f1
delta_f3dB=2;

theta0_f1=2*pi*f1/fs;
delta_theta_3dB_f1=2*pi*delta_f3dB/fs;
delta_f1=delta_theta_3dB_f1/2;
r=1-delta_f1;
c_a1_f1=r*2*cos(theta0_f1);
c_a2_f1=-r^2;
c_b0_f1=(1-r)*2*sin(theta0_f1); 

% Plotting notch filter related to carrier f1 
figure(2);
freqz([c_b0_f1],[1 -c_a1_f1 -c_a2_f1]);
title('Notch filter [carrier f1]'); % 
y_carrierf1 = filter([c_b0_f1],[1 -c_a1_f1 -c_a2_f1], y);

% Extract carrier f2
theta0_f2=2*pi*f2/fs;
delta_theta_3dB_f2=2*pi*delta_f3dB/fs;
delta_f2=delta_theta_3dB_f2/2;
r=1-delta_f2;
c_a1_f2=r*2*cos(theta0_f2);
c_a2_f2=-r^2;
c_b0_f2=(1-r)*2*sin(theta0_f2); 

% Plotting notch filter related to carrier f2 
figure(3);
freqz([c_b0_f2],[1 -c_a1_f2 -c_a2_f2]);
title('Notch filter [carrier f2]');
y_carrierf2 = filter([c_b0_f2],[1 -c_a1_f2 -c_a2_f2], y);

%% Demodulation____________________________________________________________
% Demodulation of carrier f1
y_dem1=y.*y_carrierf1;

% Demodulation of carrier f2
y_dem2=y.*y_carrierf2;

%% Bandpass filter design__________________________________________________
% Bandpass filter for carrier1
[n,fo,ao,w] = firpmord([8000 8500],[1 0],[0.0001 0.001],fs);
b = firpm(n,fo,ao,w);

% Plotting Low-Pass FIR filter frequency response 
figure(4);
freqz(b,1);
title('Low-pass FIR filter');

y_dem_bandpass1=filtfilt(b,1,y_dem1);

% Bandpass filter for carrier2
y_dem_bandpass2=filtfilt(b,1,y_dem2);

% High-pass notch IIR filter for the origin
f0_notch=0;
deltaf3dB_notch=2;
theta0_notch=f0_notch/fs*2*pi;
delta_theta3dB_notch=deltaf3dB_notch/fs*2*pi;
delta_notch=delta_theta3dB_notch/2;
r_notch=1-delta_notch;
b1_notch=-2*cos(theta0_notch);
b2_notch=-2*cos(theta0_notch);
a1_notch=2*r_notch*cos(theta0_notch);
a2_notch=-r_notch^2;

% Plotting High-pass notch IIR filter frequency response
figure(5);
freqz([1 b1_notch 1],[1 -a1_notch -a2_notch]);
title('High-pass notch IIR filter');

y_final1 = filter([1 b1_notch 1],[1 -a1_notch -a2_notch],y_dem_bandpass1);
y_final2 = filter([1 b2_notch 1],[1 -a1_notch -a2_notch],y_dem_bandpass2);

%% Demodulated signals spectrum ___________________________________________
% Spectrum of the demodulated and filtered signal #1
Y = 1/N*fftshift(fft(y_final1(:,1),N));
F = -fs/2:fs/N:fs/2-fs/N;

% Plotting spectrum of the demodulated and filtered signal
figure(6);
plot(F,abs(Y));
grid
xlabel('frequency[Hz]');
ylabel('amplitude');
title('Spectrum of the demodulated and filtered signal 1');

% Spectrum of the demodulated and filtered signal #2
Y = 1/N*fftshift(fft(y_final2(:,1),N));
F = -fs/2:fs/N:fs/2-fs/N;

% Plotting spectrum of the demodulated and filtered signal
figure(7);
plot(F,abs(Y));
grid
xlabel('frequency[Hz]');
ylabel('amplitude');
title('Spectrum of the demodulated and filtered signal 2');

%% Saving__________________________________________________________________
y_final1=y_final1.*10; % Amplifying the signal to get it louder
y_final2=y_final2.*10; % Amplifying the signal to get it louder
yfinal=[y_final1(:),y_final2(:)];
audiowrite('clean_output.wav',yfinal,fs);

