% Longitudinal force FX
% this function remap the sclar function to its vectorial form
function [fx_vec] = MF96_FX_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  fx_vec = zeros(size(kappa_vec));
  for i = 1:length(kappa_vec)

   % main code
    fx_vec(i) = MF96_FX(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
 
  end
  
 end
