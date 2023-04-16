% Lateral/longitudinal force FXY0
function [fxy0] = MF96_FXY0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [kappa__xy, Bxy, Cxy, Dxy, Exy, SVxy] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  fxy0 = magic_formula(kappa__xy, Bxy, Cxy, Dxy, Exy, SVxy);
  
 end
