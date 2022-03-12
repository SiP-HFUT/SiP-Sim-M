%%% By Wang Wenkang 2022, HFUT, China;
%%% modifed by Cheng Rui, 2022, HFUT, China

%% Descriptions:
% This function returns through- and drop-port responses of cascaded
% micro-ring resonators.

% Here it is assumed that
%   1. all the MRRs have the same radii.
%   2. the kappas (k1) between the straight bus waveguide and the ring are
%   the same.
%   3. kappas (k2)  between all the adjacent MRRs are the same.


% Variable explanations:
% g: the number of the rings;
% k1: kappa between the straight bus waveguide and ring
% k2: kappa between adjacent rings
% radius: radii of the rings

% the  parameters below should give flattop responses (for testing purpose):
%   1. g = 2 , k1=0.405, k2=0.2
%   2. g = 3 , k1 = 0.405, k2 = 0.07

function [thrs, drs] = cas_mrr (g, k1, k2, radius, betas)

l = 2*pi*radius;

% total number of the simulation wavelength points
nw = length(betas);

%% Solving symbolic expressions of the responses
syms a0 b0 a3
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

%% loop-over all the wavelength points for calculating the responses
drs = zeros(1,nw); thrs = drs;
pls = exp(-1j*betas*l/2);
T1 = tm_dc2(k1,1);
T3 = tm_dc2(k2,1);

for y = 1:nw
    pl = pls(y);
    T2 = [0 pl;pl^-1 0]; % single ring
    T4 = [0 pl;pl^-1 0]; % double ring
    Tc = T3*T2;
    Td = T1*T4*T3*T4;   
    if g == 1
        T23 = tm_dc2(k1, 1);
        T01 = T23;
        for y = 1:nw
            pl = pls(y);
            T12 = [0 pl;pl^-1 0];
            T_total = T01*T12*T23;
            [drs(y), thrs(y)] = resolve(T_total,dr_char,thr_char);
        end     
    else                          % g > 1
        if mod(g,2) == 0          % g is even
            m = g/2;
            T_total = Td^(m)*T1;          
        else                     % g is odd
            m = g;
            T_total = T1 * T2 * Tc^(m-1) * T1;
        end
        [drs(y), thrs(y)] = resolve(T_total,dr_char,thr_char);
    end
end

% DRs = 10*log10(abs(drs).^2);
% Thrs = 10*log10(abs(thrs).^2);
% figure,plot(lams*1e9,DRs),title('Drop-port response (dB)');
% figure,plot(lams*1e9,Thrs),title('through-port response (dB)');

% P_total = (abs(thrs).^2)+(abs(drs).^2);
% figure,plot(wavelengths*1e9,P_total),title('total-energy');
end

function [dr, thr] = resolve(M,dr_char,thr_char)
M1_1 = M(1,1); M1_2 = M(1,2);
M2_1 = M(2,1); M2_2 = M(2,2);
eval (['dr = ' dr_char ';'])
eval (['thr = ' thr_char ';'])
end