% Pure lateral force FY0
function [fy0] = MF96_FY0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [kappa__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  fy0 = magic_formula(kappa__y, By, Cy, Dy, Ey, SVy);
  
 end
