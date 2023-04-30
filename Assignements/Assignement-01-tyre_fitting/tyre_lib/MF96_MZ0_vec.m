<<<<<<< HEAD
 % Pure lateral force FY0
=======
% Pure aligning moment MZ0
>>>>>>> origin/main
% this function remap the sclar function to its vectorial form
function [mz0_vec] = MF96_MZ0_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  mz0_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
<<<<<<< HEAD
   % main
    mz0_vec(i) = MF96_MZ0(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
  end
  
=======
   % precode
   % main code
    mz0_vec(i) = MF96_MZ0(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

  end
  
 end
>>>>>>> origin/main
