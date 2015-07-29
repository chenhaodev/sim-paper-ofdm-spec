%%%%% 循环谱检测OFDM信号 %%%%% 
clear;
clc;
close all;

%%% OFDM参数 %%%
%数据长度的1/8必须大于等于循环谱采样长度
%数据速率 6-BPSK, 12-QPSK
%%%%%%%%%%%%%%%%
TXVECTOR.LENGTH = 256; % 数据长度
TXVECTOR.DATARATE = 6; % 数据速率
%trst_rate = 20e6; % 信号发射速率,恒定
trst_rate = 1; % 信号发射速率 normlized
%fc = 100e6;
fc = 2; % normlized

%%% 循环谱检测参数 %%%
%检测带宽为 -fs/2 至 fs/2
%循环频率分辨率为 fs/N
%频率分辨率为 M*fs/N
%%%%%%%%%%%%%%%%%%%%%
%fs = 300e6; % 采样频率
fs = 4; % 采样频率 normlized
N = 2048; % 采样长度
M = 20; % 平滑点数

%%% 信道参数 %%%
SNR = 15; % 信噪比

%{ PSK
%%% 随机数据生成 %%%
PSDU = round(rand(1,8*TXVECTOR.LENGTH));
%%% OFDM信号生成 %%%
sig = transmitter(PSDU,TXVECTOR);
ref_sig = sig;
%}

data  = randi([0 3],N,1);
sig = qammod(data,4)/sqrt(2); 
sig = sig';

%{
%%% 提升信号采样率 %%%
s_n = ceil(fs/trst_rate); % 检测采样率近似为OFDM信号速率的整数倍
sig = sig(ones(s_n,1),:); % 码元复制
sig = reshape(sig, 1, s_n*length(sig));
%}
sig = interp(sig, round(fs/fc)); %upsample

%%% 载波调制 %%%
sig_chnl = real(sig.*exp(sqrt(-1)*2*pi*fc/fs*(0:length(sig)-1)));


%%% 高斯信道 %%%
sig_awgn = awgn(sig_chnl, SNR);
%sig_awgn = awgn(zeros(1,N), SNR);

%%% 循环谱检测 %%%
cyclic_spectrum(sig_awgn, N, fs, M);
