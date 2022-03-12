% This funtion calculates the drop  and through port transmission coefficient of an admrr;
% admrr: add-drop microring resonator, whic has two bus waveguides.

function [thrs, drs] = admrr (k0, k1, radius, betas)
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

l = 2*radius*pi;

% Loop over all the wavelengths, and calculated the DR and THR responses of
% the MRR
nw = length(betas);
pls = exp(-1j*betas*l/2);
drs = zeros(1,nw); thrs = drs;
M = zeros(2,2,nw);
for i = 1:nw 
    T23 = tm_dc2 (k0,1);
    pl = pls(i);
    T12 = [0 pl; pl^-1 0];
    T10 = tm_dc2 (k1,1);
    M(:,:,i) = T10*T12*T23;
    [drs(i),thrs(i)] = resolve(M(:,:,i),dr_char,thr_char);
end

end

function [dr,thr] = resolve(M,dr_char,thr_char)
M1_1 = M(1,1); M1_2 = M(1,2);
M2_1 = M(2,1); M2_2 = M(2,2);
eval(['dr=' dr_char ';']);
eval(['thr=' thr_char ';']);
end