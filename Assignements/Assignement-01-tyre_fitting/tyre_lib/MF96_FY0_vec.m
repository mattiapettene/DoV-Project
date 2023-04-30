 % Pure lateral force FY0
% this function remap the sclar function to its vectorial form
<<<<<<< HEAD
function [fy0_vec] = MF96_FY0_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  fy0_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   % precode
   [alpha__y, By, Cy, Dy, Ey, SVy,] = MF96_FY0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
=======
function [fy0_vec, Kya_vec] = MF96_FY0_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  Kya_vec = zeros(size(alpha_vec));
  fy0_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   % precode
   [alpha__y, By, Cy, Dy, Ey, SVy, Kya, ~, ~] = MF96_FY0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
>>>>>>> origin/main
   % main code
    Kya_vec(i) = Kya;
    fy0_vec(i) = magic_formula(alpha__y, By, Cy, Dy, Ey, SVy);
  end
  