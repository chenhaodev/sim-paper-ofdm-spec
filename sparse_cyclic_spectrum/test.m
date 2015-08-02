% This example illustrates the detection for cyclostationary signal 
% band ~ [.15 .25] 
% cyclic frequency alpha = .00125.
% author: chenhaomails@gmail.com
% 2015.7

% header

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

sig.type = 'fsk'; % 'fsk'
sig.fs = 1;
sig.M = 1;

if strcmpi(sig.type,'fsk') % default signal
	load fsk.mat             
else
	error('signal type not exist!!');
end

sig.x = fsk_real;
sig.N = length(sig.x);


% test

[Spec_t, T, Tau] = cyclic_xcorr(sig.x, sig.N, sig.fs, 'show'); 

[Spec, f, alpha] = cyclic_spectrum(sig.x, sig.N, sig.fs, sig.M, 'show', 'signal_len');

