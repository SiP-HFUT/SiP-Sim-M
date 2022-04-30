clc
clearvars -except 
close all



% - Define the simulation wavelength range and the dispersion charateristic
% of the waveguide
nw = 300; lam0 = 1.55e-6; lam_span = 20e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);
alpha = 0; % lossless
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);
lams_nm = lams *1e9;

% - calculate Bragg grating response;
N=300; % total period number of the grating
k_max = 2.5e4; % max kappa of the grating
qn=ones(1,N) * k_max; % define the kappa distribution over the grating length.
period = 316e-9; 
l_g =period*(N-1); % grating length
[detuning]=lam2det(wg_info.ng,lam0,lams);
[rg,tg]=tmmcalc(qn,period,detuning);

rgi = 10*log10(abs(rg).^2);
% figure,plot(lams_nm,rgi), title('Reflection intensity response (dB) of the grating');
figure,plot(lams_nm,abs(rg)), title('Reflection coefficient (a.u) of the grating');


% - tau and kappa of the DC and MMI
% -- right MMI
t2 = sqrt(1/2);  k2 = sqrt(1/2);
% -- left DC
t1 = sqrt(1/2);  k1 = sqrt(1-t1^2);

for ps = linspace(0, pi/2, 4) 
psc = exp(-1j*ps);  
o13 = t1 * psc .* tg * sqrt(1/2) + -1j*k1 .* tg * sqrt(1/2);
o12 = t1 * psc .* rg * psc *  -1j*k1 + -1j*k1 .* rg * psc * t1;
o11 = t1 * psc .* rg * psc * t1 + -1j*k1 .* rg * psc .*  -1j*k1;

o13i = abs(o13).^2;
o12i = abs(o12).^2;
o11i = abs(o11).^2;

o_total = o11i + o12i + o13i;
figure,plot(lams_nm,o13i);
figure(22),plot(lams_nm,o_total),hold on,
end




