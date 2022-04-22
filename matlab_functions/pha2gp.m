% This function convert a phase response to a group delay (GP) response;
% based on the fact that GP response is the first derivative of the phase response

% Input Parameters: 
% w: lambdas, in m;
% pha: phase response, in rad
% w_c: center_wavelength, in m;

function [gp] = pha2gp(w,pha,lam_center)
disp('pha2gp function (lambdas (m), phase, center_wavelength (m));');
if lam_center>1
    lam_center=lam_center*1e-9;
end
if w(1)<1
    w = w*1e+9;
end
if length(pha(:,1))>1
    pha = pha';
end
ws = abs(w(2)-w(1));
delta_f = 3e+8/(lam_center)-3e+8/(lam_center+1e-9);
delta_w=ws*delta_f*2*pi/1e+12;
gp = diff(unwrap(pha))/delta_w;
gp = [gp(1) gp];
figure,plot(w,gp),title('Group delay response (ps)'),hold on;
end