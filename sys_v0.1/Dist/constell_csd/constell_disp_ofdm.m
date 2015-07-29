% This code demo the autocorrelated constellation of cyclostationary signal (e.g. ofdm) 
% author: chenhaomails@gmail.com
% 2015.7

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

sig.type = 'qam'; % 'default', 'rand', 'ofdm'
sig.snr_dB = 10; % SNR in decibels
sig.cyclic = 5; 
sig.fs = 4; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 20;
snr = 10.^(sig.snr_dB./10); % Linear Value of SNR

load signal_ofdm.mat
[s1, s2, s3] = size(signal_ofdm);
x = signal_ofdm(:,1, randi([1,100], 1,1)); x = x';
%x = (randn(1,s1)+sqrt(-1)*randn(1,s1))./sqrt(2);

r = 1; % normlized sweeping speed
M = 20; % window-size
N = length(x);

win = 'hamming';

d_tau = r/N; % time-bias resolution
tau = 0:d_tau:N; % cyclic resolution

% time-signal input 
X = x(1:N);
X = X';
i = 1; 

%% Loop
for alfa = tau

    interval_t_N = round(alfa/d_tau);
    
    % window generate
    t = 1:N;

    % spectral correlation
    X1 = X(t);
    X2 = (circshift(X(t)',interval_t_N))';
    St = conj(X1).*X2;
    St = mean(St, 1); % T average
    %S(i, floor((f_len-f_N)/2)+(1:f_N)) = St/N; %move St to central
    S(:,i) = St/N;
    i = i+1;
    
end

Spec = abs(S);

% figure

    mesh(t, tau, Spec); 
    axis tight;
    xlabel('t'); ylabel('tau');    




% autocorrelated constellation analysis
%{

cs.iter = 50;
cs.sparse = 20;
cs.ratio = 10;
cs.N = length(x);
cs.M = round(cs.N/cs.ratio);

% signal mixing
pn = pn_gen(cs.N);
mixed_sig = x.*pn;  
% ideal pass filter
lpf_temp = zeros(1,(cs.N-1)); % bandpass: side-band
lpf_temp(1) = 1;
lpf_z = interpft(lpf_temp,cs.N)/cs.N; % impulse response
% low pass filtering and downsampling
f_lpf = fft(lpf_z);
f_sig = fft(mixed_sig);
filter_out = ifft(f_lpf.*f_sig); %multiplication in freq, = convolution in time
cs.y = downsample(filter_out, cs.ratio); %defactor = downsample ratio
figure; plot(xcorr(cs.y),'o'); title('autocorrelated constellation @analog CS-RX(ratio=50), rayleigh chan, 64QAM');


%feature analysis (acf_constell)
% 1. extract/set apart the significant point from autocorrelated constellation
% 2. calculate the mean location of the left points (non-significant point,mean_g)
% 3. feature = distance(mean_g, significant point)
o1 = mean(real(cs.y));
o2 = mean(imag(cs.y));
for mm = 1:1:length(cs.y)
    temp = cs.y(mm); 
    a1 = real(temp);
    a2 = imag(temp);
    dist_point(mm) = sqrt((a1-o1)^2 + (a2-o2)^2);
end
[dist_val,dis_inx] = sort(dist_point,'descend');
key_inx = dis_inx(1);  
group_inx = setxor(dis_inx, dis_inx(1));
mean_g_r = mean(real(cs.y(group_inx))); % non-significant point
mean_g_i = mean(imag(cs.y(group_inx)));
key_g_r = real(cs.y(key_inx)); % significant point
key_g_i = imag(cs.y(key_inx));
acf_constell = sqrt((key_g_r-mean_g_r)^2 + (key_g_i-mean_g_i)^2)  %feature

%}

%{
f_mid = round(length(f)/2);
f_index_range = (f_mid - floor(length(f)/(sig.fs))+1) : (f_mid + floor(length(f)/(sig.fs)));
f_area = f(f_index_range);
a_index_range = 1 : (round(length(alpha) * 0.5));
a_area = alpha(a_index_range);
% feature extraction
y = (Spec(a_index_range,f_index_range));

figure;
mesh(f_area, a_area, Spec(a_index_range,f_index_range)); 
axis tight;
xlabel('f'); ylabel('a');
title('cyclic spectral');
norm(Spec(a_index_range,f_index_range))
%}
