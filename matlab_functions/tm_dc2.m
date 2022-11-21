%%% This function returns TYPE 2 matrix of a 2*2 directional coupler (DC)
% please refer to tm_dc2.pdf for details of this matrix

function tm_dc = tm_dc2 (kappa, tau)
tm_dc = -1/(1j*kappa).*[-tau 1; -kappa^2-tau^2 tau];
end