% Pure self-aligning moment MZ0
function [mz0] = MF96_MZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  fy0 = MF96_FY0(kappa, alpha, phi, Fz, tyre_data);
  [Br, Bt, Ct, Dr, Dt, Et, alpha__r, alpha__t, alpha__t__eq, alpha__r__eq, s] = MF96_MZ_coeffs(kappa, alpha, phi, Fz, tyre_data);
  t = pneumatic_trail(alpha__t, alpha, Bt, Ct, Dt, Et);
  mzr = residual_mzr(alpha__r, alpha, Br, Dr);

 % main code

  mz0 = -fy0 * t + mzr;
  
 end
