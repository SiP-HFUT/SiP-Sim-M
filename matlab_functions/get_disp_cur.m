% This function generates wavelength points and the corresponding effective
% indices and propagation constants at each wavelengths, from the waveguide
% and wavelength Info..

% Input Prameters
% - lam_info: struct with three fields of 'lam_center', 'lam_span' and 'nw',
%   which mean, respectivelty, center wavelength, wavelength span, and total number of the wavelength points

% - wg_info: struct with three fields of 'neff0', 'ng' and 'alpha',
%   which mean, respectivelty, effective indexand group index at the
%   wavelength of 'lam_center', wavelength span, and  light amplitude propagation loss (1/m)

% Output Results
% - lams (1*n vector): wavelength points

% - neffs (1*n vector): effective indices at each wavelength point

% - betas (1*n vector): propagation constants at each wavelength point

% Note 
% -All the wavelength units in this function are 'nm'

function [lams, neffs, betas] = get_disp_cur (lam_info, wg_info)
lam_center1 = lam_info.lam_center;
lam_span = lam_info.lam_span;
nw = lam_info.nw;
neff0 = wg_info.neff0;
ng = wg_info.ng;
lam_center2 = wg_info.lam_center;
alpha = wg_info.alpha;

lams = linspace(lam_center1-lam_span/2, lam_center1+lam_span/2, nw);
dw_2_dn = (neff0-ng)/lam_center2 ; % slope of the neff-versus-wavelengths
neff_lam = @(lam) neff0+dw_2_dn*(lam-lam_center2);
neffs = neff_lam (lams);
betas = 2*pi./(lams./neffs) - 1j*alpha; % propagation constant
end