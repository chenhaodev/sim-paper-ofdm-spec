% This code is to plot receiver operating characteristic curve for
% cs based detection for cyclostationary signal,
% , when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is pre-estimated.
% author: chenhaomails@gmail.com
% update: 15/08/05

clc; clear; close all

addpath('./Util/')
addpath('./Data/')

% Header 
snr_dB = 0; % SNR in decibels, NOTE !!!! %CS-CSD > 5
snr = 10.^(snr_dB./10); % Linear Value of SNR
gain.sig = sqrt(snr);
gain.noise = 1;

iter = 25; % Monte Carlo simulation

load bpsk.mat
mm = 1;
for thresh = 5:2:10 %CS-CSD_C, 0db
%for thresh = 0:5:20 %CS-CSD_B, 0db
%for thresh = 20:5:50 %CS-CSD_A, 0db
%for thresh = 10:2:30 %CSD
%for thresh = 20:5:50 %ED
    thresh 
    for kk=1:iter % Number of Monte Carlo Simulations
        % signal exist 
        x = bpsk(1:64); N = length(x);
        x= gain.sig.*x + gain.noise.*(randn(1,N)); 
        % signal non-exist 
        n= gain.noise.*(randn(1,N));
        
        %{
        % test with signal exist
        ed_sig(kk) = (1/length(N)).*norm(x);
        % test with signal non-exist (noise)
        ed_nos(kk) = (1/length(N)).*norm(n);
        %}
        
        %{
        % test with signal exist
        [hat_m_M1, Sx_M1, feature_mask] = disp_sparse_cyclic_spec(x);
        out1 = feature_mask.*abs(Sx_M1);
        csd_sig(kk) = (1/length(N)).*norm(out1);
        % test with signal non-exist (noise)
        [hat_m_M2, Sx_M2, feature_mask] = disp_sparse_cyclic_spec(n);
        out2 = feature_mask.*abs(Sx_M2);
        csd_nos(kk) = (1/length(N)).*norm(out2);
        %}
        
        % test with signal exist
        [hat_m_M1, Sx_M1, feature_mask] = disp_sparse_cyclic_spec(x);
        out1 = feature_mask.*abs(hat_m_M1);
        cs_csd_sig(kk) = (1/length(N)).*norm(out1);
        % test with signal non-exist (noise)
        [hat_m_M2, Sx_M2, feature_mask] = disp_sparse_cyclic_spec(n);
        out2 = feature_mask.*abs(hat_m_M2);
        cs_csd_nos(kk) = (1/length(N)).*norm(out2);
        
    end
    % Pd, Pfa
    
    %{
    ed_pd_num = length(find(ed_sig >= thresh));
    ed_pfa_num = length(find(ed_nos >= thresh));
    Pd_ed(mm) = ed_pd_num/iter;
    Pfa_ed(mm) = ed_pfa_num/iter;
    %}
    
    %{
    csd_pd_num = length(find(csd_sig >= thresh));
    csd_pfa_num = length(find(csd_nos >= thresh));
    Pd_csd(mm) = csd_pd_num/iter;
    Pfa_csd(mm) = csd_pfa_num/iter;
    %}
    
    cs_csd_pd_num = length(find(cs_csd_sig >= thresh));
    cs_csd_pfa_num = length(find(cs_csd_nos >= thresh));
    Pd_cs_csd(mm) = cs_csd_pd_num/iter;
    Pfa_cs_csd(mm) = cs_csd_pfa_num/iter;
    
    mm = mm +1;
end 
%[~, inx_ed] = sort(Pfa_ed);
%[~, inx_csd] = sort(Pfa_csd);
[~, inx_cs_csd] = sort(Pfa_cs_csd);

figure; 
%plot(sort(Pfa_ed), Pd_ed(inx_ed), '*r'); hold on;
%plot(sort(Pfa_csd), Pd_csd(inx_csd), '*g');
plot(sort(Pfa_cs_csd), Pd_cs_csd(inx_cs_csd), '*g');

%save roc_ed_0db.mat Pfa_ed Pd_ed
%save roc_csd_0db.mat Pfa_csd Pd_csd
save roc_cs_csd_0db_C.mat Pfa_cs_csd Pd_cs_csd

xlabel('Pfa');
ylabel('Pd');
