% Coefficients for Magic Formula pure longitudinal force
function [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  FZ0             = tyre_data.FZ0;
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
  LCX             = tyre_data.LCX;
  LEX             = tyre_data.LEX;
  LFZ0            = tyre_data.LFZ0;
  LGAMMAY         = tyre_data.LGAMMAY;
  LHX             = tyre_data.LHX;
  LKXK            = tyre_data.LKXK;
  LMUX            = tyre_data.LMUX;
  LVX             = tyre_data.LVX;
  

 % main code

  FZ01 = (LFZ0 * FZ0);
  dfz = Fz / FZ01 - 1;
  SHx = (dfz * pHx2 + pHx1) * LHX;
  SVx = Fz * (dfz * pVx2 + pVx1) * LVX * LMUX;
  gamma__s = (phi * LGAMMAY);
  kappa__x = kappa + SHx;
  Cx = pCx1 * LCX;
  mu__x = (dfz * pDx2 + pDx1) * (-pDx3 * gamma__s ^ 2 + 1) * LMUX;
  Dx = mu__x * Fz;
  Kxk = Fz * (dfz * pKx2 + pKx1) * exp(-(pKx3 * dfz)) * LKXK;
  Ex = (dfz ^ 2 * pEx3 + dfz * pEx2 + pEx1) * (1 - pEx4 * Sign(kappa__x)) * LEX;
  Bx = Kxk / Cx / Dx;
  
 end
