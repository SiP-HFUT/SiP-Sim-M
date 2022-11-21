function tm_dc3 = tm_dc3 (kappa, tau)
k = kappa;
t = tau;
tm_dc3 = [t 0 1j * k 0; 0 t 0 -1j * k; ... 
    1j * k 0 t 0; 0 -1j * k 0 t];
end