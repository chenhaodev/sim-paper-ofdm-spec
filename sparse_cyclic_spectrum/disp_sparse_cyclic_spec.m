% This example illustrates the cs based detection for cyclostationary signal 
% author: chenhaomails@gmail.com
% update: 15/08/05

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

% Header 

sig.type = 'fsk'; % 'fsk'
sig.fs = 1;
sig.M = 1;

if strcmpi(sig.type,'fsk') % default signal
	load fsk.mat             
else
	error('signal type not exist!!');
end

sig.x = fsk_real(1:64);
%sig.x = fsk_real(65:128);
sig.x = sig.x ./ norm(sig.x);
sig.N = length(sig.x);
n=(randn(1,sig.N)); 

% test, cs_cyc_spec + feature extract
[hat_spec] = sparse_cyclic_spec(sig.x, sig.N, sig.fs, 'show');
[out] = feature_extract(abs(hat_spec), 1:sig.N, 0.2, 1:sig.N, 0.2);
norm(out)
