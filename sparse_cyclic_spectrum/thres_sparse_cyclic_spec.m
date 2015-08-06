% calculate the threhsold in cs based detection for cyclostationary signal 
% 1. assume only noise is received, i.e., primary user is absent.
% false alarm, then
% 2. probability of false alarm = energy above threshold/No. of iteration.
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
sig.N = length(sig.x);

Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
iter = 10; % Monte Carlo simulation

%% Loop %%
for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise
		n= randn(1,sig.N); 
		sig.x = n;
		% test, cs_cyc_spec + feature extract
		[hat_spec] = sparse_cyclic_spec(sig.x, sig.N, sig.fs, 'non-show');
		[out] = feature_extract(abs(hat_spec), 1:sig.N, 0.2, 1:sig.N, 0.2);
		% energy 
        energy_fin(kk) = (1/length(sig.N)).*norm(out);% test statistic (normlized)
	end
	energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
	thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on
thres_sparse_cyclic_spec_est = thresh;
save ./Data/thres_sparse_cyclic_spec.mat thres_sparse_cyclic_spec_est

