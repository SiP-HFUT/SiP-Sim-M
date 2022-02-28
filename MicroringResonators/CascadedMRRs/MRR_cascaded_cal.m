%%% By Wang Wenkang 2022, HFUT, China;
%%% modifed by Cheng Rui, 2022, HFUT, China

%% Descriptions:
% This code returns the through- and drop-port responses of cascaded
% micro-ring resonators.

% Here it is assumed that 
%   1. all the MRRs have the same radii.
%   2. the kappas (k1) between the straight bus waveguide and the ring are
%   the same.
%   3. kappas (k2)  between all the adjacent MRRs are the same.


% variable explanations: 
    % g: the number of the rings;
    % k1: kappa between the straight bus waveguide and ring
    % k2: kappa between adjacent rings
    % radius: radii of the rings

% the parameters below should give flattop responses (for testing purpose):
%   1. g = 2 , k1=0.405, k2=0.2
%   2. g = 3 , k1 = 0.405, k2 = 0.07

clc;
close all;
clearvars except-

radius = 5e-6;

% % number of the cascaded rings
% x = input('请输入微环个数 g：');
% disp(['g = ' num2str(x)]);
% g = x;
% 
% disp('请分别输入直波导与微环之间和微环之间的kappa，k1和k2的值(中间用空格分开)：');
% k_all= str2num(input('k:','s'));
% if  g == 1
%     k1 = k_all(1);
% else
%     k1 = k_all(1);
%     k2 = k_all(2);
% end

g = 3;
% k1 = sqrt(0.45);  k2 = sqrt(0.09);
k1 = 0.405;  k2 = 0.07;

%% wavelength range and dispersion curve
sim_lam_info = struct('lam_center', 1.5546e-6, 'lam_span', 5e-9, ...
    'nw', 3000);

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2);

[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);

%% waveguide loss 
loss_per_cm_dB = 10; % typical silicon 220*500 nm waveguide loss：5 dB/cm;
loss_per_cm = 10^(loss_per_cm_dB/10);
loss_per_cm = sqrt(loss_per_cm);
alpha = log(loss_per_cm)/0.01; 
alpha = 0; % lossless;

[thrs, drs] = cas_mrr (g, k1, k2, radius, betas, alpha);

DRs = 10*log10(abs(drs).^2);
Thrs = 10*log10(abs(thrs).^2);
figure,plot(lams*1e9,DRs),title('Drop-port response (dB)');
figure,plot(lams*1e9,Thrs),title('through-port response (dB)');

