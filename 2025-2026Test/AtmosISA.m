function rho = AtmosISA(h_m, T_K)
% AtmosISA  — density (kg/m^3) using ISA pressure (with lapse rate) and ambient T
g = 9.80665; R = 287.053; T0 = 288.15; p0 = 101325; L = -0.0065;
Tstd = T0 + L*h_m;
p = p0 * (Tstd/T0)^(-g/(L*R));
rho = p / (R * T_K);
end
