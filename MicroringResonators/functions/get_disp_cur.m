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