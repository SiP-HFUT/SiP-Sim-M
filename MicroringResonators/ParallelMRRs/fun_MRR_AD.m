function M_t = fun_MRR_AD (k0,k1,radius, lams, betas,alpha)
% This funtion return the matrix, M_t, of a MRR, 
% where M_t is defined as [a3; a0] = M_t [b3; b0];

% Calculate taus of the DCs according to the assigned kappss
t0 = sqrt(1-k0^2); % upper DC
t1 = sqrt(1-k1^2); % lower DC

syms a0 b0 a3 
M =  sym('M', [2,2]);
x3 = [0;a3];
x0 = [a0;b0];
eqs = M*x3 == x0;
S = solve(eqs,[a3 b0]);
dr = simplify(S.a3/a0);
thr = simplify(S.b0/a0);
dr_char = char(dr); 
thr_char = char(thr); 

% Define the radius of the ring
l = 2*radius*pi;


% Loop over all the wavelengths, and calculated the DR and THR responses of
% the MRR
nw = length(lams);
pls = exp(-(1j*betas+alpha)*l);
drs = zeros(1,nw); thrs = drs;
M = zeros(2,2,nw);
for i = 1:nw 
    T23 = fun_TM_DC(k0,t0);
    pl = pls(i);
    T12 = [0 pl; pl^-1 0];
    T10 = fun_TM_DC(k1,t1);
    M(:,:,i) = T10*T12*T23;
    [drs(i),thrs(i)] = fun_resolve(M(:,:,i),dr_char,thr_char);
end

% DR = log10(abs(drs).^2)*10;
% THR = log10(abs(thrs).^2)*10;
% figure,plot(lams*1e9, DR),title('Drop-port response of the MRR');
% figure,plot(lams*1e9, THR),title('Through-port response of the MRR')

% M to M_t Transformation
M_t = [-M(2,1,:)./M(2,2,:)   1/M(2,2,:); 
    M(1,1,:)-M(1,2,:).*M(2,1,:)./M(2,2,:)    M(1,2,:)./M(2,2,:)];

drs1 = M_t(1,2,:)./M_t(2,2,:);
thrs1 = 1./M_t(2,2,:);
drs1 = squeeze(drs1); thrs1 = squeeze(thrs1);
DR1 = abs(drs1).^2;
THR1 = abs(thrs1).^2;
DR1_log = log10(DR1)*10;
THR1_log = log10(THR1)*10;
%figure(11),plot(lams*1e9, DR1),hold on,%title('Drop-port response of the MRR');
%figure,plot(lams*1e9, THR1),title('Through-port response of the MRR')
end
