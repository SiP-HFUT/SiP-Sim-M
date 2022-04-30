clc
clearvars -except 
close all

% - Define the simulation wavelength range and the dispersion charateristic
% of the waveguide
nw = 300; lam0 = 1.55e-6; lam_span = 15e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);
alpha = 0; % lossless
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);
lams_nm = lams *1e9;

% - calculate Bragg grating response;
k_max = 2.5e4;  % max kappa of the grating
% N=300; % total period number of the grating
% qn=ones(1,N) * k_max; % define the kappa distribution over the grating length.

N = 3000;
zn = 0:N;
k_max = 0.3e4;
qn = exp(-4*log(2)* ( (zn - ceil(N/2))/(N/2.5)).^2) * k_max; % Gaussian profile

% N = 400; % total number of the grating period
% qn=[ones(1,N/2)  -ones(1,N/2+1)] * k_max; 
figure,plot(qn);

period = 316e-9; 
l_g =period*(N-1); % grating length
[detuning]=lam2det(wg_info.ng,lam0,lams);
[rg,tg]=tmmcalc(qn,period,detuning);

rgi = 10*log10(abs(rg).^2);
% figure,plot(lams_nm,rgi), title('Reflection intensity response (dB) of the grating');
figure,plot(lams_nm,abs(rg)), title('Reflection coefficient (a.u) of the grating');


% - tau and kappa of the two MMIs
% -- right MMI
t2 = sqrt(1/2);  k2 = sqrt(1/2);
% -- left MMI
t1 = sqrt(1/2);  k1 = sqrt(1/2);

for ps = linspace(0, pi/2, 6) 
psc = exp(-1j*ps);  
o13 = t1 * tg * sqrt(1/2) + t1 .* rg * sqrt(1/2) * psc * sqrt(1/2);

o13i = abs(o13).^2;
figure (22),plot(lams_nm,o13i),hold on;

% o_total = o11i + o12i + o13i;
% figure(22),plot(lams_nm,o_total),hold on,
end




