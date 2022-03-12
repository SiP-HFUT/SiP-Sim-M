function [detuning]=lam2det(ng,lam_center,lambda)
detuning=2*pi*ng*(1./lambda-1./lam_center);
detuning=detuning(end:-1:1);
end