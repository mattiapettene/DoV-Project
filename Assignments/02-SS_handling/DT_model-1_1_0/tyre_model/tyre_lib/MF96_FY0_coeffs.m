% Coefficients for Magic Formula pure lateral force
function [alpha__y, By, Cy, Dy, Ey, SVy, Kya, SHy, mu__y] = MF96_FY0_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  Fz0             = tyre_data.Fz0;
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
  lambda__Cy      = tyre_data.lambda__Cy;
  lambda__Ey      = tyre_data.lambda__Ey;
  lambda__Fz0     = tyre_data.lambda__Fz0;
  lambda__Hy      = tyre_data.lambda__Hy;
  lambda__Kya     = tyre_data.lambda__Kya;
  lambda__Vy      = tyre_data.lambda__Vy;
  lambda__mu__y   = tyre_data.lambda__mu__y;
  lambda__gamma__y = tyre_data.lambda__gamma__y;
  

 % main code

  Fz01 = (lambda__Fz0 * Fz0);
  dfz = Fz / Fz01 - 1;
  gamma__s = (phi * lambda__gamma__y);
  SHy = (dfz * pHy2 + pHy3 * gamma__s + pHy1) * lambda__Hy;
  SVy = Fz * (pVy1 + pVy2 * dfz + (dfz * pVy4 + pVy3) * gamma__s) * lambda__Vy * lambda__mu__y;
  alpha__y = alpha + SHy;
  Cy = pCy1 * lambda__Cy;
  mu__y = (dfz * pDy2 + pDy1) * (-pDy3 * gamma__s ^ 2 + 1) * lambda__mu__y;
  Dy = mu__y * Fz;
  Ey = (dfz * pEy2 + pEy1) * (1 - (pEy4 * gamma__s + pEy3) * sign(alpha__y)) * lambda__Ey;
  Kya = Fz01 * pKy1 * sin(0.2e1 * atan((Fz / Fz01 / pKy2))) * (1 - pKy3 * abs(gamma__s)) * lambda__Fz0 * lambda__Kya;
  By = Kya / Cy / Dy;
  
 end
