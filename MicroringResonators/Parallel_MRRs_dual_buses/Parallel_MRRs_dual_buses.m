%% By Rui Cheng 2022, HFUT, China (rcheng@hfut.edu.cn)
% This code calculates drop- and through- port responses of a parallel MRR
% filter with DUAL bus waveuiges

close all;
clearvars -except p w
clc;

% define kappas (a vector) of various MRRs, from left to right; by defining
% this vector, you also have defined the number of the MRRs.
% kappas = ones(1,4)*0.2;
kappas = [0.15 0.2 0.15];


% acquire the number of the MRRs from the length of the defined kappas
nrings = length(kappas);

% radii of the MRRs, from left to right; for parallel MRRs, the MRRs generally have the same length
% radii = ones(1,nrings) * 7e-6;
radii = [7 9 7] * 1e-6;

% Define the waveguide loss
loss_per_cm_dB = 5; % typical silicon 220*500 nm waveguide loss£º5 dB/cm;
loss_per_cm = 10^(loss_per_cm_dB/10);
loss_per_cm = sqrt(loss_per_cm);
alpha = log(loss_per_cm)/0.01;
% alpha = 0; % lossless;

sim_lam_info = struct('lam_center', 1.549388e-6, 'lam_span', 0.8e-9, ...
    'nw', 2000);

% note: 'lam_center' in 'wg_info' means that 'neff0' is defined based on this
% wavelength
wg_info = struct('lam_center', 1.55e-6, 'neff0', 2.43, 'ng', 4.2);
[lams, neffs, betas] = get_disp_cur (sim_lam_info, wg_info);

lam_peak = 1.549388e-6; % observe and select the lambda of the tailored resonante peak
m = 41; % an arbitrary odd integer number

dis1 = 1.549388e-6/2.43/4*41; % distance between adjacent MRRs;
...an important parameters;
    % it should be euqual to odd number times of a quarter of the effective wavelength to ensure a resonate condition between MRRs;
% references with relevant studies need to be cited here.

for dis =[dis1]
    PL = zeros(2,2,sim_lam_info.nw);
    M = zeros(2,2,sim_lam_info.nw);
    pls = exp(-(1j*betas+alpha)*dis);
    PL(1,1,:) = pls;
    PL(2,2,:) = pls.^-1;
    
    % Define a cell that will stores the matrices of all the MRRs
    % Mi is the total transfer matrix of the i-th ring (from left to right)
    Ms = cell(1,nrings);
    
    for i = 1:nrings
        % below we will assume that the upper and lower DCs for each MRR have the
        % same kappas
        Ms{i} = tm_admrr2(kappas(i), kappas(i), radii(i), betas, alpha);
    end
    
    % Loop over all the wavelengths
    for i = 1:sim_lam_info.nw
        
        if nrings > 1
            M(:,:,i) = PL(:,:,i)*Ms{1,nrings}(:,:,i);
            for j = nrings-1 : -1 : 2
                M(:,:,i) = PL(:,:,i) * Ms{1,j}(:,:,i)* M(:,:,i);
            end
            M(:,:,i) = Ms{1,1}(:,:,i)* M(:,:,i);
        elseif nrings == 1
            M(:,:,i) = Ms{1,nrings}(:,:,i);
        end
        
    end
    drs = M(1,2,:)./M(2,2,:);
    thrs = 1./M(2,2,:);
    drs = squeeze(drs);
    thrs = squeeze(thrs);
    DR = abs(drs).^2;
    THR = abs(thrs).^2;
    DR_log = log10(abs(drs).^2)*10;
    THR_log = log10(abs(thrs).^2)*10;
    
    figure(11),plot(lams*1e9, DR),
    title(['Drop-port response of the P-MRR; disance between adjacent MRRs: ' num2str(dis*1e6) ' um']),
    hold on,
    figure(12),plot(lams*1e9, DR_log),hold on,
    figure(13),plot(lams*1e9, THR_log),hold on,
end





