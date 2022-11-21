%%% This function returns TYPE 1 matrix of a 2*2 directional coupler (DC)
% please refer to tm_dc1.pdf for details of this matrix
% t_c: transmission coefficient of the DC; for lossless DC, t_c = 1;
function tm_dc = tm_dc1 (kappa, tau)
tm_dc = [tau -1j*kappa; -1j*kappa tau];
end