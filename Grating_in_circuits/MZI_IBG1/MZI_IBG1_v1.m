% v1: modfied the code to fit the updated functions in Github (2022-03-19)
clc
clearvars -except 
close all

% - for solving the analytic expression of r and t of the system
syms a0 b0 d0 a4 c4
M =  sym('M', [4,4]);
x4 = [a4;0;c4;0];
x0 = [a0;b0;0;d0];    
eqs = M*x4 == x0; % this means we propagates the field from the output to the input
S = solve(eqs,[b0 d0 a4 c4]);
r11 = simplify(S.b0/a0);
r12 = simplify(S.d0/a0);
t13 = simplify(S.a4/a0);
t14 = simplify(S.c4/a0);
r11_char = char(r11); r12_char = char(r12); t13_char = char(t13); t14_char = char(t14);


% - Define the simulation wavelength range and the dispersion charateristic
% of the waveguide
nw = 300; lam0 = 1.55e-6; lam_span = 20e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);

% -- define loss
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


% - tau and kappa of the DCs
% -- right DC
t2 = sqrt(1/2);  k2 = sqrt(1-t2^2);
% -- left DC
t1 = sqrt(1/2);  k1 = sqrt(1-t1^2);

% --- transfer matrices of the two DCs
% TDC2 = [t1 0 1j * k1 0; 0 t1 0 -1j * k1; ... 
%     1j * k1 0 t1 0; 0 -1j * k1 0 t1];
% TDC1 = TDC2;

TDC2 = tm_dc3 (k1, 1);
TDC1 = TDC2;


labels = {};
% --- phase-shifting transfer matrix; only the upper arm has a phase shift
for ps = linspace(0, pi/2, 8) 
psc = exp(-1j*ps);    
TPS = zeros(4,4);
TPS(1,1) = psc^-1; TPS(2,2) = psc; TPS(3,3) = 1; TPS(4,4) = 1;
    
T = zeros(4, 4, nw);
r11 = zeros(1,nw); r12 = r11; t13= r11; t14 = r11;
for i = 1 : nw
% -- transfer matrices of the gratings
tgi = tg(i); rgi = 1j * rg(i);
TG = [1/tgi -rgi/tgi 0 0; rgi/tgi (tgi^2 - rgi^2)/tgi 0 0; ...
    0 0 1/tgi -rgi/tgi; 0 0 rgi/tgi (tgi^2 - rgi^2)/tgi];
    
% -- total transfer matrix
T(:,:,i) = TDC1 * TPS * TG * TDC2;

[r11(i),r12(i),t13(i),t14(i)] = ...
    resolve(T(:,:,i),r11_char,r12_char,t13_char,t14_char);
end
% figure(220),plot(lams_nm,abs(r11),'DisplayName', sprintf('ps = %.2f (pi) \n', ps/pi)),hold on,
% title(sprintf('1 to 1 ports reflection coefficient'));
% % title(sprintf('1 to 1 ports reflection response; ps = %g (pi) \n', ps/pi));


% figure(221),plot(lams_nm,abs(r12),'DisplayName', sprintf('ps = %.2f (pi) \n', ps/pi)),hold on,
% title(sprintf('1 to 2 ports reflection coefficient'));
% title(sprintf('1 to 2 ports reflection response; ps = %g (pi) \n', ps/pi));

figure,plot(lams_nm,abs(t13)),
% title(sprintf('1 to 3 ports transmission coefficient'));
title(sprintf('1 to 3 ports transmission response; ps = %g (pi) \n', ps/pi));

% o_total = abs(r11).^2 + abs(r12).^2 + abs(t13).^2 + abs(t14).^2;
% figure(221),plot(lams_nm, o_total),hold on;
end

legend ('show');

function [r11,r12,t13,t14] = resolve(M,r11_char,r12_char,t13_char,t14_char)
M1_1 = M(1,1);M1_2 = M(1,2);M1_3 = M(1,3);M1_4 = M(1,4);
M2_1 = M(2,1);M2_2 = M(2,2);M2_3 = M(2,3);M2_4 = M(2,4);
M3_1 = M(3,1);M3_2 = M(3,2);M3_3 = M(3,3);M3_4 = M(3,4);
M4_1 = M(4,1);M4_2 = M(4,2);M4_3 = M(4,3);M4_4 = M(4,4);
eval(['r11=' r11_char ';']);
eval(['r12=' r12_char ';']);
eval(['t13=' t13_char ';']);
eval(['t14=' t14_char ';']);
end