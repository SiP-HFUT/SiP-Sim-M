% - v1: modifed to use the matlab function library in github (2022-03-19)

clc
clearvars -except qn
close all

% --- for solving the analytic expression of r and t of the system
syms a0 b0 d0 a2 d2 r t
M =  sym('M', [4,4]);

x2 = [a2;r*a2+t*d2;t*a2+r*d2;d2];
% x2 = [a2;r*a2+t*d2;t*a2-r*d2;d2];
x0 = [a0;b0;0;d0];    eqs = M*x2 == x0;
S = solve(eqs,[b0 d0 a2 d2]);
r = simplify(S.b0/a0);
t = simplify(S.d0/a0);
r_char = char(r);
t_char = char(t);


% -- Define the simulation wavelength range and the dispersion
% charateristic of the waveguide, calculated from Lumerical Mode
nw = 60; lam0 = 1.55e-6; lam_span = 6e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);

% define loss
loss_db_per_cm = 5;
alpha = alpha_cal1 (loss_db_per_cm);
alpha = 0;

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);
lams_nm = lams *1e9;

% -- calculate Bragg grating response;
N=300; % define total period number of the Bragg grating
k_max = 0.95e4; % max kappa of the grating
% 
% % define the kappa distribution over the grating length.
% qn=ones(1,N) * k_max;
% % qn=[ones(1,N/2) -ones(1,N/2)] * k_max; 

N =length (qn);
period = 317e-9;

l_g =period*(N-1);
[detuning]=lam2det(wg_info.ng,lam0,lams);
[rg,tg]=tmmcalc(qn,period,detuning);

rgi = 10*log10(abs(rg).^2);
figure,plot(lams_nm,rgi), title('Reflection intensity response (dB) of the grating');
figure,plot(lams_nm,abs(rg)), title('Reflection coefficient (a.u) of the grating');


% tau and kappa of the directional coupler
t1 = sqrt(1/2);
k1 = sqrt(1-t1^2);

% length of the loop
l = 20e-6;
pls = exp(-(1j.*betas).*l/2);
re = zeros(1,nw); tran = re;tol = re;
ps = linspace(0,2*pi,16)+eps;
% ps = pi/2;
% ps = -0.75*pi;
for j = 1 : length (ps)
    
    psj = exp(-1j*ps(j));
    for i = 1:nw
        pl = pls(i);
        
        T1 = [t1 0 -1j*k1 0; 0 t1 0 1j*k1; -1j*k1 0 t1 0; 0 1j*k1 0 t1]; % left DC
       T_PS = [psj^-1 0 0 0; 0 psj 0 0; 0 0 1 0; 0 0 0 1];
%         T_PS = [1 0 0 0; 0 1 0 0; 0 0 psj^-1 0; 0 0 0 psj];
        T_L = [pl^-1 0 0 0; 0 pl 0 0; 0 0 0 pl^-1; 0 0 pl 0];
        M = T1*T_PS*T_L;
        [re(i),tran(i)] = Re(M,1j*rg(i),tg(i),r_char,t_char);
    end
    
    
    R = 10*log10((abs(re).^2));
    T = 10*log10((abs(tran).^2));
    tol = abs(re).^2+ abs(tran).^2;
figure(233),plot(lams_nm,R), title(sprintf('Reflection resposne  ps = %g (pi) \n' , ps(j)/pi)),hold on,
figure(234),plot(lams_nm,T), title(sprintf('Transmission resposne  ps = %g (pi) \n', ps(j)/pi)),hold on,
figure(335),plot(lams_nm,tol), title('Total energy'),hold on,
end
function [reflection,transmission] = Re(M,r,t,r_char,t_char)
M1_1 = M(1,1);M1_2 = M(1,2);M1_3 = M(1,3);M1_4 = M(1,4);
M2_1 = M(2,1);M2_2 = M(2,2);M2_3 = M(2,3);M2_4 = M(2,4);
M3_1 = M(3,1);M3_2 = M(3,2);M3_3 = M(3,3);M3_4 = M(3,4);
M4_1 = M(4,1);M4_2 = M(4,2);M4_3 = M(4,3);M4_4 = M(4,4);
eval(['reflection=' r_char ';']);
eval(['transmission=' t_char ';']);
end



