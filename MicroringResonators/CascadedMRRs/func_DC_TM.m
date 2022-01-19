function DC_TM = func_DC_TM(k)
    t = sqrt(1-k^2);
    DC_TM = -1/(1j*k).*[-t 1; -1 t];
end