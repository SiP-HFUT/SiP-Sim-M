% This function convert a phase response to a group delay (GP) response;
% based on the fact that GP response is the first derivative of the phase response

% Input Parameters: 
% w: lambdas, in m;
% pha: phase response, in rad
% w_c: center_wavelength, in m;

function [gp] = pha2gp(w,pha,w_c)
disp('pha2gp function (lambdas (m), phase, center_wavelength (m));');
if w_c>1
    w_c=w_c*1e-9;
end
if w(1)<1
    w = w*1e+9;
end
ws = abs(w(2)-w(1));
delta_f = 3e+8/(w_c)-3e+8/(w_c+1e-9);
delta_w=ws*delta_f*2*pi/1e+12;
gp = diff(unwrap(pha))/delta_w;
%figure(33),plot(w(2:end),gp),title('Group delay response (ps)'),hold on;
end