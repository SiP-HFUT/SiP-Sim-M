% This function convert an amplitude coefficient response to the intensity
% response in dB
% coeff: amplitude coefficient
function [I_dB] = coeff_to_IdB (coeff)
I_linear = abs(coeff).^2;
I_dB= log10(I_linear)*10;
end
