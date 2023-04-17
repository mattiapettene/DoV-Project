% Pure longitudinal force FX0
function [fx0] = MF96_FX0(kappa, alpha, phi, Fz, tyre_data, road_condition)

 % precode

  [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(kappa, alpha, phi, Fz, tyre_data, road_condition);

 % main code

  fx0 = magic_formula(kappa__x, Bx, Cx, Dx, Ex, SVx);
  
 end
