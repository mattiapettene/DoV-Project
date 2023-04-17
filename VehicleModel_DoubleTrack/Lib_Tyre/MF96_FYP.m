% Combined lateral force FY_prime (without vertical shift, used for computing combined S-A moment MZ)
function [fy_prime] = MF96_FYP(kappa, alpha, phi, Fz, tyre_data, road_condition)

 % precode

  fy0 = MF96_FY0(kappa, alpha, phi, Fz, tyre_data, road_condition);
  [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data, road_condition);

 % main code

  fy_prime = fy0 * Gyk;
  
 end
