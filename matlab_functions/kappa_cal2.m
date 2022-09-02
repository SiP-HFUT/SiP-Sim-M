%%
% delta_lambda : the bandwidths of the reflected spectrum
% ng : the group index  
% lambda_center : the central wavelength of the spectrum
% L : the length of the grating 
%%
function [kappa]=kappa_cal2(delta_lambda,ng,lambda_center,L)
kappa = sqrt((delta_lambda*pi*ng/(lambda_center)^2)^2-(pi/L)^2);
end