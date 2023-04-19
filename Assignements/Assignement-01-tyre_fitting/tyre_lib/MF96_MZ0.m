% Pure self aligning moment MZ0
function [mz0] = MF96_MZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [Br, Bt, Ct, Dr, Dt, Et, alpha__r, alpha__t] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);
                 mzr = Dr*cos(arctan(Br*alpha__r))*cos(alpha);
                 fy = MF96_FY(kappa,alpha,phi,Fz, tyre_data);
                 t = Dt*cos(Ct*arctan(Bt*alpha__t-Et(Bt*alpha__t-arctan(Bt*alpha__t))))*cos(alpha);

 % main code

  mz0 = -t * fy + mzr;
  
 end
