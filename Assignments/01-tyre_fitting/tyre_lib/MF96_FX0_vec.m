% Pure longitudinal force FX0
% this function remap the sclar function to its vectorial form
function [fx0_vec,Kxk_vec] = MF96_FX0_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  fx0_vec = zeros(size(kappa_vec));
  Kxk_vec = zeros(size(kappa_vec));
  
  for i = 1:length(kappa_vec)
   % precode
   [kappa__x, Bx, Cx, Dx, Ex, SVx, Kxk, ~, ~] = MF96_FX0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   % main code
    Kxk_vec(i) = Kxk;
    fx0_vec(i) = magic_formula(kappa__x, Bx, Cx, Dx, Ex, SVx);
  end
  
 end
