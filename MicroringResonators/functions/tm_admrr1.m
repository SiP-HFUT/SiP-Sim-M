%%% This funtion calculates Type 1 transfer matrix (illustrated in
%%% 'tm_admrr1.jpg') of an add-drop miroring resonator;

% The returned matrix has a size of 2*2*nw, where nw is the number of
% lambda points.

% admrr: add-drop microring resonator, which has two bus waveguides.
% k0,k1: kappas of the two directional couplers of the admrr
% alpha: loss per length (1/m)

function M = tm_admrr1 (k0, k1, radius, betas)

t0 = sqrt(1-k0^2); % upper DC
t1 = sqrt(1-k1^2); % lower DC

l = 2*radius*pi;

% Loop over all the wavelengths, and calculated the DR and THR responses of
% the MRR
nw = length(betas);
pls = exp(-1j*betas*l/2);
drs = zeros(1,nw); thrs = drs;
M = zeros(2,2,nw);
for i = 1:nw 
    T23 = tm_dc2 (k1, 1);
    pl = pls(i);
    T12 = [0 pl; pl^-1 0];
    T10 = tm_dc2 (k0, 1);
    M(:,:,i) = T10*T12*T23;
end

end

