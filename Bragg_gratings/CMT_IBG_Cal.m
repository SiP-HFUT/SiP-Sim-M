%  ---- Calculate the ideal reflection response of a Bragg grating according
%  to the grating strength (kappa) and phase distribution along the grating length, using
%  coupled mode theory-based TMM (CMT-TMM) method

clc
clearvars -except
close all

period = 328e-9; % grating period

% -- define the grating strength (kappa) and phase distribution along the grating length

% select one grating type and comment the others;

% - 1. uniform grating, i.e., a grating with uniform grating strength
% - along the length
% N = 300; % total number of the grating period
% qn = [ones(1,N+1)] * 4e4; % grating strength and phase distribution;  

% - 2. Gaussian-apodized grating, i.e., a grating whose grating strength
% - exhibits a Gaussian profile along the grating length
N = 800;
zn = 0:N;
qn = exp(-4*log(2)* ( (zn - ceil(N/2))/(270)).^2) * 4e4; % Gaussian profile

% - 3. pi phase-shifted grating, i.e., a grating whose grating phase
% - has an abrupt pi shift at the center of the grating
% N = 400; % total number of the grating period
% qn=[ones(1,N/2)  -ones(1,N/2+1)] * 4e4; 

figure,plot(period*(0:N)*1e3,abs(qn)/1000),title('Grating strength (*1000 m^{-1}) along the length (um)');
figure,plot(period*(0:N)*1e3,phase(qn)/pi),title('Grating phase (pi) along the length (um)');

L =period*(N);

% -- Define the simulation wavelength range and the dispersion charateristic of the waveguide, calculated from Lumerical Mode
nw = 800; lam0 = 1.55e-6; lam_span = 40e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);

% define loss
loss_db_per_cm = 5;
alpha = alpha_cal1 (loss_db_per_cm); 

wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);
lams_nm = lams *1e9;

[detuning]=lam2det(wg_info.ng,lam0,lams);
[r,t]=tmmcalc(qn,period,detuning);

ri=10*log10(abs(r).^2);
figure,plot(lams*1e+9, ri),title('Reflection intensity response (dB)')
