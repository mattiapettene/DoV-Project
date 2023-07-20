% Coefficients for Magic Formula pure longitudinal force
function [kappa__x, Bx, Cx, Dx, Ex, SVx, Kxk] = MF96_FX0_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  Fz0             = tyre_data.Fz0;
  pCx1            = tyre_data.pCx1;
  pDx1            = tyre_data.pDx1;
  pDx2            = tyre_data.pDx2;
  pDx3            = tyre_data.pDx3;
  pEx1            = tyre_data.pEx1;
  pEx2            = tyre_data.pEx2;
  pEx3            = tyre_data.pEx3;
  pEx4            = tyre_data.pEx4;
  pHx1            = tyre_data.pHx1;
  pHx2            = tyre_data.pHx2;
  pKx1            = tyre_data.pKx1;
  pKx2            = tyre_data.pKx2;
  pKx3            = tyre_data.pKx3;
  pVx1            = tyre_data.pVx1;
  pVx2            = tyre_data.pVx2;
  lambda__Cx      = tyre_data.lambda__Cx;
  lambda__Ex      = tyre_data.lambda__Ex;
  lambda__Fz0     = tyre_data.lambda__Fz0;
  lambda__Hx      = tyre_data.lambda__Hx;
  lambda__Kxk     = tyre_data.lambda__Kxk;
  lambda__Vx      = tyre_data.lambda__Vx;
  lambda__mu__x   = tyre_data.lambda__mu__x;
  lambda__gamma__y = tyre_data.lambda__gamma__y;
  

 % main code

  Fz01 = (lambda__Fz0 * Fz0);
  dfz = Fz / Fz01 - 1;
  SHx = (dfz * pHx2 + pHx1) * lambda__Hx;
  SVx = Fz * (dfz * pVx2 + pVx1) * lambda__Vx * lambda__mu__x;
  gamma__s = (phi * lambda__gamma__y);
  kappa__x = kappa + SHx;
  Cx = pCx1 * lambda__Cx;
  mu__x = (dfz * pDx2 + pDx1) * (-pDx3 * gamma__s ^ 2 + 1) * lambda__mu__x;
  Dx = mu__x * Fz;
  Kxk = Fz * (dfz * pKx2 + pKx1) * exp(-(pKx3 * dfz)) * lambda__Kxk;
  Ex = (dfz ^ 2 * pEx3 + dfz * pEx2 + pEx1) * (1 - pEx4 * sign(kappa__x)) * lambda__Ex;
  Bx = Kxk / Cx / Dx;
  
 end
