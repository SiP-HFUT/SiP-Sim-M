% This function converts wavelengths to detuning values, primarily serving as
% a preprocessing step for CMT-TMM (Coupled Mode Theory - Transfer Matrix Method)
% calculations.

function [detuning]=lam2det(ng,lam_center,lambda)
detuning=2*pi*ng*(1./lambda-1./lam_center);
detuning=detuning(end:-1:1);
end