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

sig.x = fsk_real(1:256);
sig.N = length(sig.x);
% test1, xcorr -> cyclic_spec -> cs_cyc_spec

%[Spec_cs, f, a] = cyclic_xcorr(sig.x, sig.N, sig.fs, 'show'); 

% test2, cs_cyc_spec + feature extract

[Spec_cs, f, a] = cs_cyclic_spectrum(sig.x, sig.N, sig.fs, 'show'); 
[out] = feature_test(Spec_cs, f, [-0.05 +0.05], a, [+0.01 +0.15]);
norm(out)
