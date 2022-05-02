% This function calculates kappa of a Bragg grating from its maximum
% reflectivity (0 - 1) and the length (unit: m);

% Usually, we can calculate kappa from a very short Bragg grating (as its simulation time is short).
% Then, the calculated kappa can be used in CMT-TMM for estimating 
% the spectrum features of the same Bragg grating except having a much longer length.  

% Example usage: for calculating kappa of a Bragg grating 
% with a length of 300e-6 m and a maximum reflectivity of 0.1, we can use:  
% "kappa = kappa_cal1(0.1,300e-6)"; which will give the kappa of about
% 1.09e3 m^(-1)

function [kappa]=kappa_cal1(R,L)
kappa = atanh(sqrt(R))/L;
end