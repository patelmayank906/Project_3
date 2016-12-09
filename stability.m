function dT = stability(gamma,Beta,zeta,fs)
omega = (zeta*(gamma-0.5)+sqrt((gamma/2) - Beta + zeta^2*(gamma-0.5)^2))/((gamma/2)-Beta);

dT = omega/fs(1);
end






