%% By Rui Cheng 2022, HFUT, China (rcheng@hfut.edu.cn)
% This code calculates drop- and through- port responses of a parallel MRR
% filter with a SINGLE bus waveguide

close all;
clearvars -except p w gp1
clc;

% define kappas (a vector) of various MRRs, from left to right; by defining
% this vector, you also have defined the number of the MRRs.
% kappas = ones(1,3)*0.2;
kappas = [0.32 0.32];


% acquire the number of the MRRs from the length of the defined kappas
nrings = length(kappas);

% radii of the MRRs, from left to right; for parallel MRRs, the MRRs generally have the same length
% radii = ones(1,nrings) * 7e-6;
radii = [10.4 10.4] * 1e-6;

% Define the waveguide loss
loss_per_cm_dB = 5; % typical silicon 220*500 nm waveguide loss£º5 dB/cm;
alpha = alpha_cal1(loss_per_cm_dB);
% alpha = 0; % lossless;

sim_lam_info = struct('lam_center', 1.5455e-6, 'lam_span', 40e-9, ...
    'nw', 5000);

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);



dis = 10e-6; % distance between adjacent MRRs;
...unlike parallel MRRs with dual-bus, here dis WON'T impact the thransmission response of the system.
    
pls = exp(-1j*betas*dis);

% Define a cell, t_mrrs, that will stores the transmission coefficient responses of all the MRRs
% ti is the transmission coefficient response of the i-th ring (from left to right)
t_mrrs = cell(1,nrings);

for dw = linspace(0.0,1,5) % dw: wavelength shift
pss = [dw -dw];
for i = 1:nrings
    t_mrrs{i} = t_mrr_ap(kappas(i), radii(i), betas, pss(i));
end

t = t_mrrs{1};
for j = 2 : 1 : nrings
    t = pls .* t_mrrs{j} .* t;
end

ithr = log10(abs(t).^2)*10;
pthr = phase(t);
gp = pha2gp(lams,pthr,1.55e-6);
figure,plot(lams*1e9, ithr),
title('Through-port intensity response of the P-MRR');
figure,plot(lams(2:end)*1e9, gp),
title('Through-port group delay response (ps) of the P-MRR');
end






