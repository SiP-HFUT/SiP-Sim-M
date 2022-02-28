%%% This function returns the transmission coefficient versus lambda, t(lambdas), of an all-pass MRR (i.e., MRR with a single bus waveguide),

% kap: kappa of the directional coupler;
% radius = radius of the ring;
% lams: a scalar listing all the lambdas of the simulation;
% betas: a scalar listing all the propgation constants; = 2pi*neff/lambda;
% alpha: propagation loss ( /m )

function thrs = t_mrr_ap (kap, radius, betas, alpha)

tau = sqrt(1-kap^2); 

syms a0 b0 a1 b1 pml
M =  sym('M', [2,2]);
x1 = [a1;b1];
x1 = subs (x1, b1, a1 * pml);
x0 = [a0;b0];
eqs = M*x1 == x0;
S = solve(eqs,[a1 b0]);
thr = simplify(S.b0/a0);
thr_char = char(thr);

l = 2*radius*pi;

% Loop over all the wavelengths, and calculated the DR and THR responses of
% the MRR
nw = length(betas);
pls = exp(-(1j*betas+alpha)*l);
thrs = zeros(1,nw);
for i = 1:nw
    T = tm_dc2 (kap, 1);
    M = T;
    thrs(i) = resolve(M, pls(i)^-1, thr_char);
end

% for validation purpose
% i_thr = log10(abs(thrs).^2)*10;
% pha_thr = phase(thrs) / pi;
% figure,plot(lams*1e9, i_thr),title('Through-port intensity response of the MRR');
% figure,plot(lams*1e9, pha_thr),title('Through-port phase response of the MRR');
end

function [thr] = resolve(M, pml, thr_char)
M1_1 = M(1,1); M1_2 = M(1,2);
M2_1 = M(2,1); M2_2 = M(2,2);
eval(['thr=' thr_char ';']);
end
