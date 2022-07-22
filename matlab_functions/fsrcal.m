% This function calculates FSR of a microring resonator
% dl (unit: m): for MZIs, it is the length difference between the two arms; for MRRs, it is the perimeter of the ring  
% ng: group index of the waveguide
% lam0 (unit: m): center wavelength of the interested wavelength band
function [fsr] = fsrcal (dl, ng, lam0)
fsr = lam0^2 / (ng*dl);
end
