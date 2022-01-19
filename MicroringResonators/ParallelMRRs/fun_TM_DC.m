% function to define the transfer matrix of a DC; TM: transfer matrix;
function TM_DC = fun_TM_DC(k,t)
    TM_DC = -1/(1j*k).*[-t 1; -1 t];
end