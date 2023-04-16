% Self aligning moment
function [mz0] = MF96_FMZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [kappa__mz, Bmz, Cmz, Dmz, Emz, SVmz] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  mz0 = magic_formula(kappa__m, Bm, Cm, Dm, Em, SVm);
  
 end
