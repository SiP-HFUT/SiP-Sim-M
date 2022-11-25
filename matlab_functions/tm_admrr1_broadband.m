%%% This function is the same as tm_admrr1.m, except that it is used for
%%% broadband, i.e., alpha, kappa and tau are wavelength-dependent, which means that each of these parameters is a 1*nw matrix (nw is the number of the wavelength points), rather than single
% values

function M = tm_admrr1_braodband (k0s,t0s, k1s,t1s, radius, betas, ps)

l = 2*radius*pi;

% Loop over all the wavelengths, and calculated the DR and THR responses of
% the MRR
nw = length(betas);
pls = exp(-1j*betas*l/2) .* exp(-1j*ps/2);
drs = zeros(1,nw); thrs = drs;
M = zeros(2,2,nw);
for i = 1:nw 
    T23 = tm_dc2 (k1s(i), t1s(i));
    pl = pls(i);
    T12 = [0 pl; pl^-1 0];
    T10 = tm_dc2 (k0s(i), t0s(i));
    M(:,:,i) = T10*T12*T23;
end

end

