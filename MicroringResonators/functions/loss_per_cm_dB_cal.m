function loss_per_cm_dB = loss_per_cm_dB_cal(alpha)
loss_per_cm_dB = log10((exp(alpha*0.01))^2)*10;
fprintf('The loss per cm is about %g dB \n', loss_per_cm_dB);
end