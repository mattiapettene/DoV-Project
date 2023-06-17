% Combined self-aligning moment MZ0
function [mz] = MF96_MZ(kappa, alpha, phi, Fz, tyre_data)

 % precode

  fy_prime = MF96_FYP(kappa, alpha, phi, Fz, tyre_data);
  fx = MF96_FX(kappa, alpha, phi, Fz, tyre_data);
  [Br, Bt, Ct, Dr, Dt, Et, alpha__r, alpha__t, alpha__t__eq, alpha__r__eq, s] = MF96_MZ_coeffs(kappa, alpha, phi, Fz, tyre_data);
  t__eq = pneumatic_trail(alpha__t__eq, alpha, Bt, Ct, Dt, Et);
  mzr__eq = residual_mzr(alpha__r__eq, alpha, Br, Dr);

 % main code

  mz = s * fx - fy_prime * t__eq + mzr__eq;
  
 end
