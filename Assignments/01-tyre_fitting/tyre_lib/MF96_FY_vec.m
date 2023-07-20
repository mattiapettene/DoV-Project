% Longitudinal force FX
% this function remap the sclar function to its vectorial form
function [fy_vec,Gyk_vec] = MF96_FY_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  fy_vec = zeros(size(kappa_vec));
  Gyk_vec = zeros(size(kappa_vec));
  for i = 1:length(kappa_vec)

   % main code
    [fy_vec(i),Gyk_vec(i)] = MF96_FY(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
 
  end
  
 end
