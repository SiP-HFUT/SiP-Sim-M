%% By Rui Cheng 2022, HFUT, China (rcheng@hfut.edu.cn)
% This code calculate the through and drop responses of an add-drop
% microring resonator (admrr)

% Difference with ad_mrr.m: this code uses wavelength-dependent alpha, kappa
% and tau, which means that each of these parameters is a 1*nw matrix (nw is the number of the wavelength points), rather than single
% values

clc; % clear the command window
close all; % close all the figures plotted before
clearvars -except 

%% First, Ddefine the simulation wavelength range and the dispersion charateristic of the waveguide, calculated from Lumerical Mode
nw = 10000; % first define the total number of the wavelength points.
os = ones (1,nw); % only for convenience
sim_lam_info = struct('lam_center', 1.55e-6, 'lam_span', 60e-9, ...
    'nw', nw);

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);



%% Define the parameters of the ring
radius = 10e-6; % radius of the ring
l = 2*radius*pi;

% Define the waveguide loss
% - direct definition
a = 1; % a: transmission coefficient per round trip of the ring
alpha = alpha_cal2 (a, l); 
alphas = os .* alpha;
% - imported from external data and then perform interpolation
alphas = interp1 (lamsnm_ex, alphas_ex, lams*1e9);

% Define the DC parameters
% - direct definition
% t0s = os * 0.98 * 0.98; t1s = os * 0.98;
% k0s = sqrt(1-t0s.^2); k1s = sqrt(1-t1s.^2); 
% k0s = os * 0.1; k1s = os * 0.1;
% t0s = sqrt(1-k0s.^2); t1s = sqrt(1-k1s.^2); 

% - imported from external data and then perform interpolation
t0s = interp1 (lamsnm_ex, t0s_ex, lams*1e9);
k0s = interp1 (lamsnm_ex, k0s_ex, lams*1e9);
t1s = interp1 (lamsnm_ex, t1s_ex, lams*1e9);
k1s = interp1 (lamsnm_ex, k1s_ex, lams*1e9);


%% for solving the analytic expression of drop- and through-port responses of the system
syms a0 b0 a3
M =  sym('M', [2,2]);
x3 = [0;a3];
x0 = [a0;b0];
eqs = M*x3 == x0;
S = solve(eqs,[a3 b0]);
dr = simplify(S.a3/a0);
thr = simplify(S.b0/a0);
dr_char = char(dr);
thr_char = char(thr);





% set a phase shift inside the ring to shift the resonance wavelength if neccesary
phase_shift = 0; 
M = tm_admrr1_broadband (k0s, t0s, k1s, t1s, radius, betas, phase_shift);
drs = zeros(1,sim_lam_info.nw); thrs = drs;
for i = 1 : sim_lam_info.nw
    [drs(i),thrs(i)] = resolve(M(:,:,i), dr_char, thr_char);
end

DR = log10(abs(drs).^2)*10;
THR = log10(abs(thrs).^2)*10;
figure(1),plot(lams*1e9, DR),title('Drop-port response of the MRR');hold on;
figure(2),plot(lams*1e9, THR),title('Through-port response of the MRR');hold on;

pthrs = phase(thrs);
figure,plot(lams*1e9, pthrs/pi),title('Through-port phase response (pi)');

% for validation purpose
p_total = abs(drs).^2 + abs(thrs).^2;
figure,plot(lams*1e9, p_total),title('Total energy');

function [dr,thr] = resolve(M,dr_char,thr_char)
M1_1 = M(1,1); M1_2 = M(1,2);
M2_1 = M(2,1); M2_2 = M(2,2);
eval(['dr=' dr_char ';']);
eval(['thr=' thr_char ';']);
end

