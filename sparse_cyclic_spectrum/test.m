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
%sig.x=(randn(1,sig.N)+sqrt(-1)*randn(1,sig.N))./(sqrt(2)); % complex noise

% test2, cs_cyc_spec + feature extract

[Spec_cs, Spec, f, a] = cs_cyclic_spectrum(sig.x, sig.N, sig.fs, 'show'); 
[out1] = feature_test(Spec, f, [-0.2 +0.2], a, [-0.2 +0.2]);
[out2] = feature_test(Spec_cs, f, [-0.2 +0.2], a, [-0.2 +0.2]);
norm(out1)
norm(out2)
