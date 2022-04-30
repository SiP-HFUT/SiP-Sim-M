clc
clearvars -except 
close all

%% for solving the analytic expression of r and t of the system
%if exist ('r_char') ==0
syms a0 b0 d0 a2 d2 r t
M =  sym('M', [4,4]);

x2 = [a2;r*a2+t*d2;t*a2+r*d2;d2];
x0 = [a0;b0;0;d0];    eqs = M*x2 == x0;
S = solve(eqs,[b0 d0 a2 d2]);
r = simplify(S.b0/a0);
t = simplify(S.d0/a0);
r_char = char(r);
t_char = char(t);                                                                                                                                                                                                                                                                                                                                   
%end
%return



%% calcualte the response of the IBG (PhC) reflector;
span =20e-9;
nw=600
wavelength0=1550e-9;
wavelengths = wavelength0 + linspace(-0.5*span, 0.5*span, nw);
ng = 4.2;
wnm = wavelengths*1e+9;

N=300;
%zn=1:N;
%qn=exp(-4*log(2)* ( (zn - N/2 )/(270) ).^2); % Gaussian apodized
qn=[ones(1,N)];
period = 316e-9;
L_IBG =period*(N-1);


% %Ideal spectrum
[detuning]=lam2det(ng,wavelength0,wavelengths); % ORIGIANL N_G=4.195;
[r_g,t_g]=tmmcalc(qn*9e+3,period,detuning);   % max(abs(qnG))/max(abs(qn))
R_G = 10*log10(abs(r_g).^2);
T_G = 10*log10(abs(t_g).^2);
figure,plot(wnm,R_G),title('Reflection of the reflector (IBG or PhC)');
figure,plot(wnm,T_G),title('Transmission of the reflector (IBG or PhC)');

%return

radius = 10e-6;
 l = 2*pi*radius;%l2 = 10e-6;l3 = 10e-6;
% a = 1;
% alpha = log(0.5)/-0.01;
% exp(-alpha*l1)

w0 = 1.55;
neff0 = 2.43;
dw_2_dn = (neff0-ng)/w0;
neff_wavelength = @(w) neff0+dw_2_dn*(w*1e+6-w0);
neff_wavelength = neff_wavelength (wavelengths);
belta = 2*pi./(wavelengths./neff_wavelength);

for  delta_phi = 0:pi/6:pi/2
    ps = exp(-1j*delta_phi);
  



%% Calculate MRR response without the gratings
t1 = sqrt(1/2);
k1 = sqrt(1-t1^2);

re = zeros(1,nw); tran = re;tol = re;
for i = 1:nw
    pl = exp(-(1j*belta(i)).*l);
    T1 = [t1 0 -1j*k1 0; 0 t1 0 1j*k1; -1j*k1 0 t1 0; 0 1j*k1 0 t1]; % left DC
    T_PS = inv([ps 0 0 0; 0 ps^-1 0 0; 0 0 1 0; 0 0 0 1]);
    T_L = [pl^-1 0 0 0; 0 pl 0 0; 0 0 0 pl^-1; 0 0 pl 0];
    M = T1*T_PS*T_L;
    [re(i),tran(i)] = Re(M,1j*r_g(i),t_g(i),r_char,t_char);
end


R = 10*log10((abs(re).^2));
T = 10*log10((abs(tran).^2));
tol = abs(re).^2+ abs(tran).^2;
figure(12),plot(wnm,R), title('Reflection resposne'),hold on,
figure(13),plot(wnm,T), title('Transmission resposne'),hold on,
figure(14),plot(wnm,tol), title('Total energy'),hold on,
end
function [reflection,transmission] = Re(M,r,t,r_char,t_char)
M1_1 = M(1,1);M1_2 = M(1,2);M1_3 = M(1,3);M1_4 = M(1,4);
M2_1 = M(2,1);M2_2 = M(2,2);M2_3 = M(2,3);M2_4 = M(2,4);
M3_1 = M(3,1);M3_2 = M(3,2);M3_3 = M(3,3);M3_4 = M(3,4);
M4_1 = M(4,1);M4_2 = M(4,2);M4_3 = M(4,3);M4_4 = M(4,4);
eval(['reflection=' r_char ';']);
eval(['transmission=' t_char ';']);
end



