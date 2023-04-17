% Combined longitudinal force FX
function [fx] = MF96_FX(kappa, alpha, phi, Fz, tyre_data, road_condition)

 % precode

  fx0 = MF96_FX0(kappa, alpha, phi, Fz, tyre_data, road_condition);
  [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data, road_condition);

 % main code

  fx = fx0 * Gxa;
  
 end
