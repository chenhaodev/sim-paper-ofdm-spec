% This code plots the ROC curve for Rayleigh channel, based on the paper
% titled "On the Energy Detection of Unknown Signals
% over Fading Channels." In particular, the plot corresponds to Fig. 1 of
% the paper.
% This code is written by Sanket S. Kalamkar, Indian Institute of
% Technology, Kanpur.
clc
close all
clear all
L = 10; % Number of sensing samples to be taken
snr_db = 20; % Average SNR in decibel for Rayleigh channel
snr = 10.^(snr_db./10);
thresh = 0:0.01:100; % Threhsold
% Calculation of probability of false alarm
Pf= 1- gammainc(thresh./2, L./2); % the 'upper' incomplete gamma function
pd = [];
A = snr./(1 + snr);
u = L./2; % Time-Bandwidth product
for pp = 1:length(Pf)
n = 0:1:u-2;
term_sum1 = sum((1./factorial(n)).*(thresh(pp)./2).^(n));
term_sum2 = sum((1./factorial(n)).*(((thresh(pp)./2).*(A)).^(n)));
pd(pp) = exp(-thresh(pp)./2).*term_sum1 + (1./A).^(u-1).*(exp(-thresh(pp)./(2.*(1+snr))) - exp(-thresh(pp)./2).*term_sum2); % Probability of detection % don't get the equation
end

% roc curve fig1 on ref paper
%loglog(Pf,1-pd,'r') 
%grid on
%xlim([10^-4 1])
%ylim([10^-5 1])
% roc curve normal
plot(Pf, pd, 'r')
grid on

