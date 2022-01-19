%% By Rui Cheng 2022, HFUT, China (rcheng@hfut.edu.cn)

% This code calculates the drop- and through-port responses of a Fano
% resonator based on dual MRR interference.
% This kind of Fano resonators can be found in Fig. 5.4 of Hu Ting's PhD thesis.  
% The final output of such Fano resonators is actually given by the combination of one
% MRR's drop-port light with another MRR's through-port light

close all
clearvars except-
clc

% Define the waveguide loss
loss_per_cm_dB = 5; % typical silicon 220*500 nm waveguide loss£º5 dB/cm;
loss_per_cm = 10^(loss_per_cm_dB/10);
loss_per_cm = sqrt(loss_per_cm);
alpha = log(loss_per_cm)/0.01; 
% alpha = 0; % lossless;

sim_lam_info = struct('lam_center', 1.5500374e-6, 'lam_span', 1e-9, ...
    'nw', 1000);

% note: 'lam_center' in 'wg_info' means that the constant 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2);
[lams, neffs, betas] = fun_get_disp_cur (sim_lam_info, wg_info);
lams_nm = lams*1e9;
% radii of the two MRRs; 
radius1 = 10e-6;
radius2 = radius1;

% kappas of the DC of two MRRs; 
kap1 = 0.2;
kap2 = 0.2;


%  dphi: the phase difference between the two interfered lights
for dphi = pi

% calculate the responses of the two MRRs   
[drs1, thrs1] = fun_MRR_AD(kap1, kap1, radius1, lams, betas,alpha);
[drs2, thrs2] = fun_MRR_AD(kap2, kap2, radius2, lams, betas,alpha);


%% for Output 1 (upper output in Fig. 5.4)
% expressions of the two interfered lights
light1 = sqrt(1/2) * thrs1; 
light2 = sqrt(1/2) * -1j * exp(-1j * dphi) * drs2;

figure,plot(lams*1e9, abs(light1).^2), title('Amplitude response of the light1');
figure,plot(lams_nm, phase(light1)/pi), title('Phase response of the light1 (pi)');
 
figure,plot(lams*1e9, abs(light2).^2), title('Amplitude response of light2');
figure,plot(lams_nm, phase(light2)/pi), title('Phase response of light2 (pi)');

dphase = (phase(light1) - phase(light2)) / pi;
figure,plot(lams_nm, dphase), 
title('Phase difference between Lights 1 and 2 (pi) versus lambda');
figure,plot(lams_nm, 0.5 * cos(dphase).^2), 
title('Amplitude modulation profile due to the phase difference between Lights 1 and 2');

% The out light is essentially the superposition of lights 1 and 2.
out1 = (light1 + light2); 
am_out1 = 10 * log10 (abs(out1).^2); % amplitude response
figure,plot(lams_nm, am_out1), title(['Output 1 response: dphi = ' num2str(dphi/pi) ' pi']);


%% for Output 2 (lower output in Fig. 5.4)
% % expressions of the two interfered lights
% light1 = sqrt(1/2) * drs1; 
% light2 = sqrt(1/2) * -1j * exp(-1j * dphi) * thrs2;
% out2 = (light1 + light2); 
% am_out2 = 10 * log10 (abs(out2).^2); % amplitude response
% figure,plot(lams_nm, am_out2), title(['Output 2 response; dphi = ' num2str(dphi/pi) ' pi']);

% % calculate the total energy of the two outputs
% figure, plot(lams_nm, abs(out1).^2 + abs(out2).^2 ), 
% title('Total energy of the two outputs');
end





