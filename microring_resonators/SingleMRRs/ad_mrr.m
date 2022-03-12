%% By Rui Cheng 2022, HFUT, China (rcheng@hfut.edu.cn)
% This code calculate the through and drop responses of an add-drop
% microring resonator (admrr)

clc; % clear the command window
close all; % close all the figures plotted before

%% Define the parameters of the ring
radius = 10e-6; % radius of the ring
l = 2*radius*pi;

% Define the waveguide loss
% a = 0.96; % a: transmission coefficient per round trip of the ring
% alpha = alpha_cal2 (0.96, l); 

% or define loss per cm in dB, typically 5 dB for 200*500 nm SOI waveguide
loss_db_per_cm = 5;
alpha = alpha_cal1 (loss_db_per_cm); 

t0 = 0.97; t1 = 0.9735; k0 = sqrt(1-t0^2); k1 = sqrt(1-t1^2); 
% k0 = 0.1; k1 = 0.1;

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



%% Define the simulation wavelength range and the dispersion charateristic of the waveguide, calculated from Lumerical Mode
sim_lam_info = struct('lam_center', 1.55e-6, 'lam_span', 30e-9, ...
    'nw', 6000);

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);

M = tm_admrr1 (k0, k1, radius, betas);
drs = zeros(1,sim_lam_info.nw); thrs = drs;
for i = 1 : sim_lam_info.nw
    [drs(i),thrs(i)] = resolve(M(:,:,i), dr_char, thr_char);
end

DR = log10(abs(drs).^2)*10;
THR = log10(abs(thrs).^2)*10;
figure,plot(lams*1e9, DR),title('Drop-port response of the MRR');
figure,plot(lams*1e9, THR),title('Through-port response of the MRR');

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

