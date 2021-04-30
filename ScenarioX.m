clearvars -except Hd Hd1
close all
clc

[data_orig, fs] = audioread('original.wav'); %original signal in t
data_fft_orig = fft(data_orig); %original signal in f
[data_scr, fs] = audioread('scrambled.wav'); %scrambled signal in t
data_fft_scr = fft(data_scr); %scrambled signal in f

fs = 44100;            % Sampling frequency                    
ts = 1/fs;             % Sampling period
N=length(data_orig);   % Number of samples 
%t = (1:N)*ts;         % Time vector
f = horzcat(-linspace(0,N/2,N/2)*fs/N,linspace(N/2,0,N/2)*fs/N);

%Here is the original signal
subplot(1,2,1)
plot(f,abs(data_fft_orig))
title('Original Signal in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 9e3])
ylim([0 800])

%Here is the scrambled signal
subplot(1,2,2)
plot(f,abs(data_fft_scr))
title('Scrambled Signal in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 9e3])
ylim([0 800])

%Alright. We have made a lowpass filter called "Hd1". Let's implement it:
data_filtered=filter(Hd1,data_scr);
data_fft_filtered=fft(data_filtered);

%Here is the scrambled signal after filtering
figure
plot(f,abs(data_fft_filtered))
title('Filtered Signal in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 9e3])
ylim([0 800])

%Here is the impulse response of the filter
figure
[h,t]=impz(Hd1,fs); %t=44100 or 44099 from this.
t=t*1000/fs; %adjust time to actually be time.
plot(t,h) %comparing with filterDesigner, exactly the same.
title('Impulse Response of Filter')
xlabel('Time (ms)')
ylabel('Amplitude')
xlim([0 30])
ylim([-0.13 0.25])

%Here is the frequency response of the filter.
%***I can't seem to derive it from the impulse response.***
figure
[H,w]=freqz(Hd1,fs);
w=w*fs/1000/2/pi; %adjust frequency to actually be frequency
plot(w,20*log10(abs(H))) %comparing with filterDesigner, exactly the same.
title('Frequency Response of Filter')
xlabel('Frequency (kHz)')
ylabel('Gain (dB)')
xlim([0 max(w)])
ylim([-105 1])

% H=fft(h); %this kind of works
% plot(w,20*log10(H))

fs = 44100; % Sampling frequency                    
ts = 1/fs; % Sampling period
N=220500; %220500 samples
t = ts:ts:5; % Time vector
f = horzcat(-linspace(0,N/2,N/2)*fs/N,linspace(N/2,0,N/2)*fs/N);

figure
sig=sin(2*pi*2000*t)+sin(2*pi*4000*t)+sin(2*pi*6000*t)+sin(2*pi*8000*t)+sin(2*pi*10000*t)+sin(2*pi*12000*t);
plot(t(1:200),sig(1:200));
title('Sample Signal in Time Domain')
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 200*ts])
ylim([-6 6])

figure
sig_trans=fft(sig);
plot(f,abs(sig_trans))
xlim([0 1.4e4])
title('Sample Signal in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

figure
filt_sig=filter(Hd1,sig);
filt_sig_trans=fft(filt_sig);
plot(f,abs(filt_sig_trans))
xlim([0 1.4e4])
title('Filtered Sample Signal in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

figure
plot(t(1:200),filt_sig(1:200));
title('Filtered Sample Signal in Time Domain')
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 200*ts])
ylim([-6 6])
%%
%failures
% bb = ones(1, L) / L;
% ww = -pi:(pi / 100):pi;
% HH = freqz(bb, 1, ww);
% 
% subplot(2, 1, 1), plot(ww, abs(HH));
% subplot(2, 1, 2), plot(ww, angle(HH));
% figure
% H=0; %initializing H
% z=tf('z',ts); %
% for n=1:length(h)
%     H=H+h(n)*z^-n;
% end
%The process
% [data_orig, fs] = audioread('original.wav'); %original signal in t
% data_fft_orig = fft(data_orig); %original signal in f
% %My code
% % subplot(1,2,1)
% plot(f,abs(data_fft_orig))
% title('Original Signal in Frequency Domain')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% xlim([0 1e4])
% ylim([0 1200])
% hold on
% % %Matlab code
% % P2 = abs(data_fft_orig/L); %two-sided spectrum P2
% % P1 = P2(1:L/2+1); %one-sided spectrum P1
% % P1(2:end-1) = 2*P1(2:end-1); 
% % f = fs*(0:(L/2))/L; %frequency domain f
% % subplot(2,2,3)
% % plot(f,P1) 
% % title('Single-Sided Amplitude Spectrum of X(t)')
% % xlabel('f (Hz)')
% % ylabel('|P1(f)|')
% 
% [data_scr, fs] = audioread('scrambled.wav'); %scrambled signal in t
% data_fft_scr = fft(data_scr); %scrambled signal in f
% %My code
% % subplot(1,2,2)
% plot(f,abs(data_fft_scr))
% title('Scrambled Signal in Frequency Domain')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% xlim([0 1e4])
% ylim([0 1200])

% %Matlab code
% P2_scr = abs(data_fft_scr/L); %two-sided spectrum of scrambled signal
% P1_scr = P2_scr(1:L/2+1); %one-sided spectrum of scrambled signal
% P1_scr(2:end-1) = 2*P1_scr(2:end-1);
% f = fs*(0:(L/2))/L; %frequency domain of scrambled signal
% subplot(2,2,4)
% plot(f,P1_scr) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

% figure

%ok this doesn't quite work...why don't i try to find gain for each
%frequency where magnitude is greater than 10??
% plot(abs(data_fft_scr(:,1))./abs(data_fft_orig(:,1)))
% title('Gain')
% xlim([1.75e5 2.22e5])

%figure

%plot(pwelch(data_orig))

% figure
% 
% X=zeros(1,220500);
% for n=0:220499
%     X(1,1+n)=abs(data_fft_scr(220500-n))./abs(data_fft_orig(1+n));
% end
% 
% x=1:220500;
% plot(x,X(1,x))
% xlim([1.75e5 2.22e5])
% title('')
% ylim([0 2.22e5])