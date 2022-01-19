%% By Rui Cheng, HFUT, China (rcheng@hfut.edu.cn)
% This code calculates the through port resposne of a 4-port optical router with WG crossing and dual
% MRRs, whose schematic is the same as that in Fig. 4.19 (Cross State) in Hu Ting's PhC thesis.  

close all;
clearvars except-
clc;

% Define the waveguide loss
loss_per_cm_dB = 0; % typical silicon 220*500 nm waveguide loss£º5 dB/cm;
loss_per_cm = 10^(loss_per_cm_dB/10);
loss_per_cm = sqrt(loss_per_cm);
alpha = log(loss_per_cm)/0.01; 
% alpha = 0; % lossless;

sim_lam_info = struct('lam_center', 1.55e-6, 'lam_span', 3e-9, ...
    'nw', 2600);

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2);
[lams, neffs, betas] = fun_get_disp_cur (sim_lam_info, wg_info);

% radii of the two MRRs; 
radius1 = 10e-6;
radius2 = radius1;

% kappas of the DC of two MRRs; 
kap1 = 0.2;
kap2 = 0.2;

% the aditional optical path(OP), OP(I1-to-O1) - OP(I1-to_R2-O1);
% here we simply assume that these two have a dphi (rad) difference in
% phase
for dphi = 1.4%0:0.2:pi/2

% calculate the transfer matrix (TM) of the bottom-left (BL) MRR   
[drs_bl, thrs_bl] = fun_MRR_AD(kap1, kap1, radius1, lams, betas,alpha);

% TM of the top-right (TR) MRR  
[drs_tr, thrs_tr] = fun_MRR_AD(kap2, kap2, radius2, lams, betas,alpha);

% output of the two paths of the light
light1 = drs_bl;
light2 = thrs_bl .* drs_tr * exp(-1j * dphi);

% figure,plot(lams*1e9, abs(drs_bl).^2), title('Amplitude response of the 1st path of light');
figure,plot(lams*1e9, phase(drs_bl)/pi), title('Phase response of the 1st path of light (pi)');

% figure,plot(lams*1e9, abs(light2).^2), title('Amplitude response of the 2rd path of light');
figure,plot(lams*1e9, phase(light2)/pi), title('Phase response of the 2rd path of light (pi)');

figure,plot(lams*1e9, (phase(drs_bl) + phase(light2)) / pi), title('Summary phase of the two path of the light (pi) versus lambda');
figure,plot(lams*1e9, cos (phase(drs_bl) + phase(light2)/2).^2); 
% The response of I1 is essentially the one superimposed by two paths of light.
o1 = (light1 * sqrt(0.5) + light2 * -1j * sqrt(0.5)); 

P_o1 = 10 * log10 (abs(o1).^2);
Phase_o1 = phase(o1)/pi;
figure,plot(lams*1e9, P_o1), title(['dphi = ' num2str(dphi) ' rad']);
%figure,plot(lams*1e9, Phase_o1), title(['dphi = ' num2str(dphi) ' rad']);

% % phase response
% pha_o1 = phase (res_o1)/pi;
% figure,plot(lams*1e9, pha_o1), title('Phase response (pi)');
% 
% % group delay response
% gp_o1 = pha2gp(lams,pha_o1,sim_lam_info.lam_center);
% figure,plot(lams(2:end) * 1e9, gp_o1),title('Group delay response (ps)');
end





