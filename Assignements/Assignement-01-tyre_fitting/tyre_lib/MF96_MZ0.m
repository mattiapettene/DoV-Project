% Pure self aligning moment MZ0
function [mz0] = MF96_MZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode
 % all arctan are substituited with atan MODIFY MAPLE!

  [Br, Bt, Ct, Dr, Dt, Et, alpha__r, alpha__t] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);
                 mzr = Dr*cos(atan(Br*alpha__r))*cos(alpha);
                 fy0 = MF96_FY0(kappa,alpha,phi,Fz, tyre_data);
                 t = Dt*cos(Ct*atan(Bt*alpha__t-Et*(Bt*alpha__t-atan(Bt*alpha__t))))*cos(alpha); %sistemare il et"*"(Bt... in maple

 % main code

  mz0 = -t * fy0 + mzr;
  
 end
