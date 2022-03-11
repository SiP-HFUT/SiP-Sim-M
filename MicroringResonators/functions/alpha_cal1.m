function alpha = alpha_cal1 (loss_per_cm_dB)
loss_per_cm = 10^(loss_per_cm_dB/10);
loss_per_cm = sqrt(loss_per_cm);
alpha = log(loss_per_cm)/0.01; 
end
