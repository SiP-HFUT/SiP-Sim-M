clc
clearvars -except 
close all

%% for solving the analytic expression of r and t of the system
syms a0 b0 d0 a3 c3
M =  sym('M', [4,4]);
x3 = [a3;0;c3;0];
x0 = [a0;b0;0;d0];    eqs = M*x0 == x3;
S = solve(eqs,[b0 d0 a3 c3]);
r11 = simplify(S.b0/a0);
r12 = simplify(S.d0/a0);
t13 = simplify(S.a3/a0);
t14 = simplify(S.c3/a0);
r11_char = char(r11);
r12_char = char(r12);
t13_char = char(t13);
t14_char = char(t14);

%% calcualte the response of the IBG (PhC) reflector;
nw = 200; lam0 = 1.55e-6; lam_span = 20e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);
% define loss
alpha = 0; % lossless
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);
lams_nm = lams *1e9;

N=200; % define total period number of the Bragg grating
k_max = 3e4; % max kappa of the grating
qn=[ones(1,N)] * k_max; % define the kappa distribution over the grating length.
period = 316e-9;
l_g =period*(N-1);
[detuning]=lam2det(wg_info.ng,lam0,lams);
[rg,tg]=tmmcalc(qn,period,detuning);

rgi = 10*log10(abs(rg).^2);
figure,plot(lams_nm,rgi), title('Reflection intensity response (dB) of the grating');
figure,plot(lams_nm,abs(rg)), title('Reflection coefficient (a.u) of the grating');

for  delta_phi = 1*pi%:pi/3:pi
    ps = exp(-1j*delta_phi);

dr = zeros(1,nw);
thr = dr;


t1 = sqrt(1/2);
k1 = sqrt(1-t1^2);
t2 = sqrt(1/2);
k2 = sqrt(1-t2^2);
r11 = zeros(1,nw); r12 = r11; t13 = r11; t14 = r11;t_tot = r11;
for i = 1:nw
    T1 = [t1 0 -1j*k1 0; 0 t1 0 1j*k1; -1j*k1 0 t1 0; 0 1j*k1 0 t1]; % left DC
    rei = re(i);tri = tr(i);
    
    %rei = abs(re(i));tri = abs(tr(i));
    %tri = eps; rei = 1-eps;  % r an t for the IBG
    
    % from [r t;r t], together with 1j*r;
    rei = 1j*rei;
    T2 = [tri-rei^2/tri rei/tri 0 0; -rei/tri 1/tri 0 0; 0 0 tri-rei^2/tri rei/tri; 0 0 -rei/tri 1/tri];
    
    % from [r t; conj(t) -conj(r)]. without 1j*r;
    %T2 = [1/tri rei/tri 0 0; rei/tri 1/tri 0 0; 0 0 1/tri rei/tri; 0 0 rei/tri 1/tri];
    
    T3 = T1;
    T_PS = [ps 0 0 0; 0 ps^-1 0 0; 0 0 1 0; 0 0 0 1];
    M = T3*T2*T_PS*T1;
    [r11(i),r12(i),t13(i),t14(i)] = Re(M,r11_char,r12_char,t13_char,t14_char);
end


R11 = 10*log10((abs(r11).^2));
R12 = 10*log10((abs(r12).^2));
T13 = 10*log10((abs(t13).^2));
T14 = 10*log10((abs(t14).^2));
figure(12),plot(wnm,R12), title('Reflection resposne 1-to-2 port'),hold on,
figure(13),plot(wnm,R11), title('Reflection resposne 1-to-1 port'),hold on,
figure(14),plot(wnm,T14), title('Transmission resposne 1-to-4 port'),hold on,
figure(15),plot(wnm,T13), title('Transmission resposne 1-to-3 port'),hold on,
%t_tot = abs(t13).^2+abs(t14).^2;figure,plot(t_tot);
end
function [r11,r12,t13,t14] = Re(M,r11_char,r12_char,t13_char,t14_char)
M1_1 = M(1,1);M1_2 = M(1,2);M1_3 = M(1,3);M1_4 = M(1,4);
M2_1 = M(2,1);M2_2 = M(2,2);M2_3 = M(2,3);M2_4 = M(2,4);
M3_1 = M(3,1);M3_2 = M(3,2);M3_3 = M(3,3);M3_4 = M(3,4);
M4_1 = M(4,1);M4_2 = M(4,2);M4_3 = M(4,3);M4_4 = M(4,4);
eval(['r11=' r11_char ';']);
eval(['r12=' r12_char ';']);
eval(['t13=' t13_char ';']);
eval(['t14=' t14_char ';']);
end



