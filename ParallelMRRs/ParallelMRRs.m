%% By Rui Cheng, HFUT, China (rcheng@hfut.edu.cn) 
% This code calculate drop- and through- port responses of a parallel MRR filter

close all;
clearvars except-
clc;

imported_funs_MRR = funs_MRR();

k0 = 0.1;
k1 =k0;
radius = 7e-6; % radius of the ring

% Define the waveguide loss
% loss_per_cm_dB = 2.5; % typical loss: 5 dB per cm;
% loss_per_cm = -log10(loss_per_cm_dB/10);
% loss_per_cm = sqrt(loss_per_cm);
% alpha = -log(loss_per_cm)/0.01; % typical waveguide loss£º6 dB/cm;
alpha = 0; % lossless;

sim_lam_info = struct('lam_center', 1.549387e-6, 'lam_span', 1e-9, ...
    'nw', 1000);
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2);
[lams, neffs, betas] = fun_get_disp_cur (sim_lam_info, wg_info); 

dis1 = 1.549387e-6/2.43/4*41; % distance between adjacent MRRs; 
...an important parameters to ensure resonate conditions 
    
for dis =[dis1]
PL = zeros(2,2,sim_lam_info.nw);
pls = exp(-(1j*betas+alpha)*dis);
PL(1,1,:) = pls;
PL(2,2,:) = pls.^-1;
M_total = PL;



M1 = imported_funs_MRR.fun_MRR_AD(k0, k1, radius, lams, betas,alpha);
M2 = M1;
M3 = M1;
for i = 1:sim_lam_info.nw
M_total(:,:,i) = M1(:,:,i)*PL(:,:,i)*M2(:,:,i)*PL(:,:,i)*M3(:,:,i);
end
drs = M_total(1,2,:)./M_total(2,2,:);
thrs = 1./M_total(2,2,:);
drs = squeeze(drs);
thrs = squeeze(thrs);
DR = abs(drs).^2;
THR = abs(thrs).^2;
DR_log = log10(abs(drs).^2)*10;
THR_log = log10(abs(thrs).^2)*10;

figure(11),plot(lams*1e9, DR),
title(['Drop-port response of the P-MRR; disance= ' num2str(dis*1e6) ' um']),
hold on,
end
%figure,plot(lams*1e9, THR),title('Through-port response of the P-MRR')





