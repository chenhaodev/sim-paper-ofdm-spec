% Here we calculate the threhsold in energy detection by simulations. 
% This is a general method and applicable to all scenarios for energy detection.
% We assume that all the signals are complex Gaussian.
% Algorithm:
% 1. Assume only noise is received, i.e., primary user is absent.
% false alarm, then
% 2. Probability of False Alarm = energy above threshold/No. of Iteration.
% Code written by: Sanket Kalamkar, Indian Institute of Technology Kanpur,
% India.

clc
close all
clear all
L = 100; % 
iter = 100; % number of Monte Carlo Simulations
Pf = 0.01:0.01:1; % Probability of False Alarm
for tt = 1:length(Pf)
    tt
	for kk=1:iter 
		% noise1
		%n=(randn(1,L)+j*randn(1,L))./(sqrt(2)); % complex, real + image
		% noise2
		n=randn(1,L); % real
		y = n; 
		% energy
		energy = abs(y).^2; 
		energy_fin(kk) =(1/L).*sum(energy); % test statistic (normlized)
	end
 energy_desc = sort(energy_fin,'descend'); % arrange values in descending order
 thresh(tt) = energy_desc(ceil(Pf(tt)*iter)); % largest thres to avoid 'Pf(tt)' in each 'iter' group
end
plot(thresh, Pf, 'ob')
hold on

%%
thresh1 = (qfuncinv(Pf)./sqrt(L))+ 1; % Theroretical value of threshold
plot(thresh1, Pf, 'r')
hold on
