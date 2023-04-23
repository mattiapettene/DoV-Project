% Pure longitudinal force FX0
% this function remap the sclar function to its vectorial form
function [fx0_vec] = MF96_FX0(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  fx0_vec = zeros(size(kappa_vec));
  for i = 1:length(kappa_vec)
   % precode
   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   % main code
    fx0_vec(i) = magic_formula(kappa__x, Bx, Cx, Dx, Ex, SVx);
  end
  
 end
