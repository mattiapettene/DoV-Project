% Combined longitudinal force FX
function [fx] = MF96_FX(kappa, alpha, phi, Fz, tyre_data)

 % precode

  fx0 = MF96_FX0(kappa, alpha, phi, Fz, tyre_data);
  [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  fx = fx0 * Gxa;
  
 end
