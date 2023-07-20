% Combined lateral force FY
function [fy] = MF96_FY(kappa, alpha, phi, Fz, tyre_data)

 % precode

  fy0 = MF96_FY0(kappa, alpha, phi, Fz, tyre_data);
  [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  fy = fy0 * Gyk + SVyk;
  
 end
