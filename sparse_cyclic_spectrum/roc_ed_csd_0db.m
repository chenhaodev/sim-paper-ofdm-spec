clc; clear; close all

addpath('./Util/')
addpath('./Data/')


load roc_ed_0db.mat
load roc_csd_0db.mat
load roc_cs_csd_0db.mat

[~, inx_ed] = sort(Pfa_ed);
[~, inx_csd] = sort(Pfa_csd);
[~, inx_cs_csd] = sort(Pfa_cs_csd);

plot(sort(Pfa_ed), Pd_ed(inx_ed), '*r'); hold on;
plot(sort(Pfa_csd), Pd_csd(inx_csd), '*g'); hold on;
plot(sort(Pfa_cs_csd), Pd_cs_csd(inx_cs_csd), '.b'); 

xlabel('Pfa');
ylabel('Pd');
legend('ED', 'CSD', 'CS-CSD');
title('ROC of ED, CSD, CS-CSD');