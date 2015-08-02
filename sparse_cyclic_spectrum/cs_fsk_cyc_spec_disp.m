% This code demo the cyclic spectrum density of cyclostationary signal (e.g. ofdm) 
% author: chenhaomails@gmail.com
% 2015.7

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

sig.type = 'fsk'; % 'default', 'rand', 'ofdm'
sig.cyclic = 16; 
sig.fs = 4; %normlized sample rate 
sig.fc = 2; %normlized RF-carrier rate
sig.M = 5;

sig.snr_dB = -10; % SNR in decibels
snr = 10.^(sig.snr_dB./10); % Linear Value of SNR

load FSK_real.mat
x = fsk_real(1:512);
%noise 
%x =(randn(1,512)+sqrt(-1)*randn(1,512))./(sqrt(2));
%signal + noise
%x = fsk_real(1:512)+(randn(1,512)+sqrt(-1)*randn(1,512))./(sqrt(2));

% cs sampling and sparse reconstruction 
cs.sparse = 16;
cs.ratio = 4;
cs.iter = 10;
cs.N = length(x);
cs.M = round(cs.N/cs.ratio);
Phi = randn(cs.M,cs.N);
y = Phi*x';
%A = Phi*(dftmtx(cs.N))^-1;
%A = Phi*(haarmtx(cs.N))^-1;
A = Phi*(dctmtx(cs.N))^-1;

[recov,~] = cosamp(y, A, cs.sparse, cs.iter);  
[max_value,max_inx] = max(abs(recov));
min_inx = find(recov < abs(max_value) * 0.01); % threshold
recov(min_inx) = 0;
%hatx = (haarmtx(cs.N))^-1 * recov;
%hatx = (dftmtx(cs.N))^-1 * recov;
hatx = (dctmtx(cs.N))^-1 * recov;


		% standard cyclic spectrum detection
		[Spec_o, f_o, alpha_o] = my_cyclic_spectrum(x, cs.N, sig.fs, sig.M,'show', 'standard_cyclic_spectrum');
		% cs based cyclic spectrum detection
		[Spec, f, alpha] = my_cyclic_spectrum(hatx', cs.N, sig.fs, sig.M,'show', 'cs_based_cyclic_spectrum');
		%feature range: f ~[-fs/3 fs/3]; a ~[0, a/2]
		f_mid = round(length(f)/2);
		f_index_range = (f_mid - floor(length(f)/2/3)+1) : (f_mid + floor(length(f)/2/3));
		f_area = f(f_index_range);
		a_index_range = 1 : (ceil(length(alpha) * 0.5));
		a_area = alpha(a_index_range);
		% feature extraction
		y = (Spec(a_index_range,f_index_range));
		y_o = (Spec_o(a_index_range,f_index_range));
		% main component
		%{
		Spec_r = reshape(Spec, 1, length(f)* length(alpha));
		sorted_spec = sort(Spec_r, 'descend'); 
		main_spec = sorted_spec(1:10);
		rest_spec = sorted_spec(11:end);
		test_stat = norm(main_spec) ./ norm(rest_spec)
		%}
