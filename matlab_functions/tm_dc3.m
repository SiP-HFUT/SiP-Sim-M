function tm_dc3 = tm_dc3 (kappa, tc)
k = tc * kappa;
t = sqrt (tc^2 - k^2);
tm_dc3 = [t 0 1j * k 0; 0 t 0 -1j * k; ... 
    1j * k 0 t 0; 0 -1j * k 0 t];
end