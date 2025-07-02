% coefficient to intensity (dB) conversion
function [y] = c2i(x)
y = log10(abs(x).^2)*10;
end