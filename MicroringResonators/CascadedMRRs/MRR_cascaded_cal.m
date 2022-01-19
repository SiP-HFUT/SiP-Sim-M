%% By Wang Wenkang 2022, HFUT, China

%% Descriptions:
% This code calculates the through- and drop-port responses of cascaded
% micro-ring resonators.
% variable explanations: 
    % g: the number of the rings;
    % k1: kappa between the straight bus waveguide and ring
    % k2: kappa between adjacent rings

%%include func_DC_TM
clc;
close all;
clearvars except-

%% ��������
syms a0 b0 a3 g k1 k2;
M = sym('M',[2,2]);  
if mod(g,2) == 0   % even
    x0 = [a0;b0];
    x3 = [0;a3];
    eqs = M*x3 == x0;
    S = solve(eqs,[a3;b0]);
    dr = S.a3/a0;
    thr = S.b0/a0;
    dr_char = char(dr);
    thr_char = char(thr);

else               % odd
    x0 = [b0;a0];
    x3 = [0;a3];
    eqs = M*x3 == x0;
    S = solve(eqs,[a3;a0]);
    dr = S.a3/b0;
    thr = S.a0/b0;
    dr_char = char(dr);
    thr_char = char(thr);

end

%% part two ΢������
radius = 10e-6;
l = 2*pi*radius;

%% part three ������沨����Χ
wavelength_center = 1549.26e-9;
span = 1e-9;
nw = 600;
wavelengths = linspace(wavelength_center-span/2,wavelength_center+span/2,nw);
%% 500*220nm   ��������
neff = 2.43;
ng = 4.2;
slope_neff = (neff-ng)/1559e-9;
neff_wavelength = @(wavelength) neff + slope_neff * (wavelength - 1559e-9 );

neff_wavelengths = neff_wavelength(wavelengths);
betas = 2*pi ./ wavelengths .* neff_wavelengths;


%% part four define waveguide loss�������
loss_per_cm_db = 0;
loss_per_cm = 10^(-loss_per_cm_db/10);
alpha = -log(sqrt(loss_per_cm))/0.01; 
% alpha = 0;

%% ����΢������
x = input('������g��');
disp(['g = ' num2str(x)]);
g = x;


%% ���岨����K��T  
% ���Բ���  ring = 2 , k1=0.405 k2=0.2
% ring = 3 , k1 = 0.405 k2 = 0.07
disp('������k1��k2��ֵ(�м��ÿո�ֿ�)��');
k_all= str2num(input('k:','s'));
if  g == 1
    k1 = k_all(1);
else
 
    k1 = k_all(1);
    k2 = k_all(2);
end

%% part six ѭ��������в���������ȡ�˿ں�ֱͨ�˿ڵ���Ӧ
drs = zeros(1,nw); thrs = drs;
pls = exp(-1j*betas*l/2)*exp(-alpha*l/2);
T1 = func_DC_TM(k1); 
T3 = func_DC_TM(k2);

for y = 1:nw
    pl = pls(y);
    T2 = [0 pl;pl^-1 0]; %single ring 
    T4 = [0 pl;pl^-1 0]; %double ring
    Tc = T3*T2;
    Td = T1*T4*T3*T4;

    if g == 1
    %% COPY FROM testmrr.m
    T23 = func_DC_TM(k1);
    T01 = T23;
    for y = 1:nw
     %   beta=betas(i);
        pl = pls(y);
        T12 = [0 pl;pl^-1 0];
        T_total = T01*T12*T23;
        [drs(y), thrs(y)] = Rsolve(T_total,dr_char,thr_char);
    end

    else                          % g > 1
        if mod(g,2) == 0          % even   gΪż��  
            m = g/2;
            T_total = Td^(m)*T1;

         else                     % odd   gΪ����
            m = g;
            T_total = T1 * T2 * Tc^(m-1) * T1;
        end
         [drs(y), thrs(y)] = Rsolve(T_total,dr_char,thr_char);
    end
end

%% 
figure,plot(abs(drs).^2);
DRs = 10*log10(abs(drs).^2);
Thrs = 10*log10(abs(thrs).^2);
figure,plot(wavelengths*1e9,DRs),title('Drop-port response');
figure,plot(wavelengths*1e9,Thrs),title('through-port response');

% P_total = (abs(thrs).^2)+(abs(drs).^2);
% figure,plot(wavelengths*1e9,P_total),title('total-energy');
%%
function [dr, thr] = Rsolve(M,dr_char,thr_char)
M1_1 = M(1,1); M1_2 = M(1,2);
M2_1 = M(2,1); M2_2 = M(2,2);
eval (['dr = ' dr_char ';'])
eval (['thr = ' thr_char ';'])
end