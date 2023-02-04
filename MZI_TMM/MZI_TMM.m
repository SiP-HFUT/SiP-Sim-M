% This codes calcualte transsmission resposne of a MZI, using transfer matrix method 
clc
clearvars -except 
close all

% - define the simulation wavelength range and the dispersion charateristic
% of the waveguide
nw = 600; lam0 = 1.55e-6; lam_span = 60e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);

% - define loss
alpha = 0; % lossless

% - dispersion characteristic

neff0 = 2.43;
wg_info = struct('lam_center', 1.55e-6, 'neff0', neff0, 'ng', 4.2, 'alpha', alpha);
[lams, neffs0, betas0] = get_disp_cur (sim_lam_info, wg_info);
lamsnm = lams *1e9;


for dt = linspace (0, 43, 5) % ig you have temperature variations
neff1 = 2.43 + dt*1.8e-4;
wg_info = struct('lam_center', 1.55e-6, 'neff0', neff1, 'ng', 4.2, 'alpha', alpha);
[lams, neffs1, betas1] = get_disp_cur (sim_lam_info, wg_info);

% - define kappas and taus of the DCs
k1 = sqrt(1/2);
t1 = sqrt (1-k1^2);
k2 = sqrt(1/2);
t2 = sqrt (1-k2^2);

% - define the length of the top and bottom arms
l1 = 100e-6;
l2 = l1 + 30e-6;

T0 = [1; 0];
T1 = tm_dc1 (k1, t1);
T3 = tm_dc1 (k2, t2);
thr = zeros(1,nw);
cro = thr;
for i = 1:nw
T2 = [exp(-1j*betas1(i)*l1) 0; 0 exp(-1j*betas0(i)*l2)];
T_total = T3 * T2 * T1 * T0;
thr(i) = T_total(1);
cro(i) = T_total(2);
end

ithr = abs(thr).^2;
icro = abs(cro).^2;

lgthr = log10 (ithr) * 10;
lgcro = log10 (icro) * 10;

figure,plot(lamsnm, lgthr),hold on, plot(lamsnm, lgcro),
title(['delta T = ' num2str(dt)]);
legend('through (dB)', 'cross (dB)')
end





