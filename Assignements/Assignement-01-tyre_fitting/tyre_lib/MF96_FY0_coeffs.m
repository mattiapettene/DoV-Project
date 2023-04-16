% Coefficients for Magic Formula pure lateral force
function [kappa__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  FZ0             = tyre_data.FZ0;
  Fz01            = tyre_data.Fz01;
  pCy1            = tyre_data.pCy1;
  pDy1            = tyre_data.pDy1;
  pDy2            = tyre_data.pDy2;
  pDy3            = tyre_data.pDy3;
  pEy1            = tyre_data.pEy1;
  pEy2            = tyre_data.pEy2;
  pEy3            = tyre_data.pEy3;
  pEy4            = tyre_data.pEy4;
  pHy1            = tyre_data.pHy1;
  pHy2            = tyre_data.pHy2;
  pHy3            = tyre_data.pHy3;
  pKy1            = tyre_data.pKy1;
  pKy2            = tyre_data.pKy2;
  pKy3            = tyre_data.pKy3;
  pVy1            = tyre_data.pVy1;
  pVy2            = tyre_data.pVy2;
  pVy3            = tyre_data.pVy3;
  pVy4            = tyre_data.pVy4;
  LCY             = tyre_data.LCY;
  LEY             = tyre_data.LEY;
  LFZ0            = tyre_data.LFZ0;
  LGAMMAY         = tyre_data.LGAMMAY;
  LHY             = tyre_data.LHY;
  LKA             = tyre_data.LKA;
  LMUY            = tyre_data.LMUY;
  LVY             = tyre_data.LVY;
  

 % main code

  FZ01 = (LFZ0 * FZ0);
  dfz = Fz / FZ01 - 1;
  gamma__s = (phi * LGAMMAY);
  SHy = (dfz * pHy2 + pHy3 * gamma__s + pHy1) * LHY;
  SVy = Fz * (pVy1 + pVy2 * dfz + (dfz * pVy4 + pVy3) * gamma__s) * LVY * LMUY;
  alpha__y = alpha + SHy;
  Cy = pCy1 * LCY;
  mu__y = (dfz * pDy2 + pDy1) * (-pDy3 * gamma__s ^ 2 + 1) * LMUY;
  Dy = mu__y * Fz;
  Ey = (dfz * pEy2 + pEy1) * (1 - (pEy4 * gamma__s + pEy3) * Sign(alpha__y)) * LEY;
  Kya = Fz01 * pKy1 * sin(0.2e1 * atan((Fz / Fz01 / pKy2))) * (1 - pKy3 * my_abs(gamma__s)) * LFZ0 * LKA;
  By = Kya / Cy / Dy;
  
 end
