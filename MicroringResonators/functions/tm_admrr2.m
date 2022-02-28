%%% This funtion calculates Type 2 transfer matrix (illustrated in
%%% 'tm_admrr2.pdf') of an add-drop miroring resonator;

% The returned matrix has a size of 2*2*nw, where nw is the number of
% the lambda points.

% admrr: add-drop microring resonator, which has two bus waveguides.
% k0,k1: kappas of the two directional couplers of the admrr
% alpha: loss per length (1/m)

function M = tm_admrr2 (k0, k1, radius, betas,alpha)

M0 = tm_admrr1 (k0, k1, radius, betas, alpha);

% Type 1 to  Type 2 Matrix Transformation
M = [-M0(2,1,:)./M0(2,2,:)   1/M0(2,2,:);
    M0(1,1,:)-M0(1,2,:).*M0(2,1,:)./M0(2,2,:)    M0(1,2,:)./M0(2,2,:)];

end
