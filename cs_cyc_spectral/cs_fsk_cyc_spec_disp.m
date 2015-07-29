% This code demo the cyclic spectrum density of cyclostationary signal (e.g. ofdm) 
% author: chenhaomails@gmail.com
% 2015.7

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

sig.type = 'fsk'; % 'default', 'rand', 'ofdm'
sig.snr_dB = 10; % SNR in decibels
sig.cyclic = 16; 
sig.fs = 4; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 5;
snr = 10.^(sig.snr_dB./10); % Linear Value of SNR

load FSK_real.mat
x = fsk_real(1:512);
%noise 
%x =(randn(1,512)+sqrt(-1)*randn(1,512))./(sqrt(2));

% cs sampling and sparse reconstruction 
cs.sparse = 8;
cs.ratio = 8;
cs.iter = 64;
cs.N = length(x);
cs.M = round(cs.N/cs.ratio);
Phi = randn(cs.M,cs.N);
y = Phi*x';
A = Phi*(dftmtx(cs.N))^-1;
[recov,~] = cosamp(y, A, cs.sparse, cs.iter);  
hatx = (dftmtx(cs.N))^-1 * recov;

% standard cyclic spectrum detection
[Spec, f, alpha] = my_cyclic_spectrum(x, cs.N, sig.fs, sig.M,'show', 'standard_cyclic_spectrum');
% cs based cyclic spectrum detection
[Spec_cs, f_cs, alpha_cs] = my_cyclic_spectrum(hatx', cs.N, sig.fs, sig.M,'show', 'cs_based_cyclic_spectrum');
% main component
Spec_r = reshape(Spec, 1, length(f_cs)* length(alpha_cs));
sorted_spec = sort(Spec_r, 'descend'); 
main_spec = sorted_spec(1:10);
rest_spec = sorted_spec(11:end);
test_stat = norm(main_spec) ./ norm(rest_spec)