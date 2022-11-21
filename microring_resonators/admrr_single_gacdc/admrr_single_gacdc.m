% ---- calculate the drop- and through-port responses of a mrr with one
% GA-CDC and a regular directional coupler as the light couplers; by Cheng Rui, 2022; 

clc
clearvars -except qn 
close all

% -- Define the simulation wavelength range and the dispersion charateristic of the waveguide, calculated from Lumerical Mode
nw = 1200; lam0 = 1.55e-6; lam_span = 30e-9;
sim_lam_info = struct('lam_center', lam0, 'lam_span', lam_span, ...
    'nw', nw);

% define loss
loss_db_per_cm = 5;
alpha = alpha_cal1 (loss_db_per_cm); 

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2, 'alpha', alpha);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);
lams_nm = lams *1e9;

% -- calculate integrated Bragg grating (IBG) response;
N=300; % define total period number of the IBG
qn=[ones(1,N)];

%zn=1:N;
%qn=exp(-4*log(2)* ( (zn - N/2 )/(270) ).^2); % Gaussian apodized

period = 316e-9;
l_g =period*(N-1);
k_max = 0.2e4; % max kappa of the grating
[detuning]=lam2det(wg_info.ng,lam0,lams); % ORIGIANL N_G=4.195;
[r,t]=tmmcalc(qn*k_max,period,detuning);   % max(abs(qnG))/max(abs(qn))

ri = 10*log10(abs(r).^2);
figure,plot(lams_nm,ri), title('Reflection intensity response (dB) of the grating')
figure,plot(lams_nm,abs(r)), title('Reflection coefficient (a.u) of the grating')


% -- define the ring parameters
% kappa of the directional coupler
tau = 0.9;  kappa=sqrt(1-tau^2);
l_total = 200e-6; % l_total: total length of the ring (including the grating);
l_ng = l_total - l_g; % l_ng: length of the non-grating part of the ring
if l_ng <=0
    fprintf('Error! the ring is not long enough to accommodate the grating! \n');
return;
end
% -- solve the analytic expressions of the drop- and through-port responses of the system
syms a0 b0 b3 
M =  sym('M', [2,2]);
x3 = [0;b3];
x0 = [a0;b0];
eqs = M*x3 == x0;
S = solve(eqs,[b3 b0]);
thr = simplify(S.b0/a0);
drop = simplify(S.b3/a0);
drop_char = char(drop);
thr_char = char(thr);

% add a phase shift to the ring to fine tune the resonate wavelength to match the Bragg resonance of the grating
for ps = [0.573]*pi 
p1 = exp(-1j*betas.*l_ng/2).*exp(-1j*ps);

dr = zeros(1,nw);
thr = dr;

% -- Calculate MRR response with the gratings (using its complex response)
for j = 1:nw
    rj = 1j*r(j); tj = t(j);
    T1 = [1/rj -tj/rj; tj/rj (rj^2-tj^2)/rj];
    T2 = [p1(j)^-1 0; 0 p1(j)];
    T3 = tm_dc2 (kappa, tau);
    M = T1*T2*T3;
    [dr(j),thr(j)] = Re(M,drop_char,thr_char);
end

THR = 10*log10((abs(thr).^2));
DR = 10*log10((abs(dr).^2));
% figure,plot(lams_nm,THR), title(sprintf('Through-port intensity response (dB); ps = %g (pi) \n', ps/pi));
figure,plot(lams_nm,DR), title(sprintf('Drop-port intensity response (dB); ps = %g (pi) \n', ps/pi));
end

function [r,t] = Re(M,drop_char,thr_char)
M1_1 = M(1,1);M1_2 = M(1,2);%M1_3 = M(1,3);M1_4 = M(1,4);
M2_1 = M(2,1);M2_2 = M(2,2);%M2_3 = M(2,3);M2_4 = M(2,4);
%M3_1 = M(3,1);M3_2 = M(3,2);M3_3 = M(3,3);M3_4 = M(3,4);
%M4_1 = M(4,1);M4_2 = M(4,2);M4_3 = M(4,3);M4_4 = M(4,4);
eval(['r=' drop_char ';']);
eval(['t=' thr_char ';']);
end



