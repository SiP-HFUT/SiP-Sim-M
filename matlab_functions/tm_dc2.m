%%% This function returns TYPE 2 matrix of a 2*2 directional coupler (DC)
% please refer to tm_dc2.pdf for details of this matrix
% t_c: transmission coefficient of the DC; for lossless DC, t_c = 1;
function tm_dc = tm_dc2 (kappa, t_c)
kappa = t_c * kappa;
tau = sqrt (t_c^2 - kappa^2);
tm_dc = -1/(1j*kappa).*[-tau 1; -t_c^2 tau];
end