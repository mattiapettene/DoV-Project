% Coefficients for Magic Formula combined forces
function [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  Fz0             = tyre_data.Fz0;
  rBx1            = tyre_data.rBx1;
  rBx2            = tyre_data.rBx2;
  rBy1            = tyre_data.rBy1;
  rBy2            = tyre_data.rBy2;
  rBy3            = tyre_data.rBy3;
  rCx1            = tyre_data.rCx1;
  rCy1            = tyre_data.rCy1;
  rHx1            = tyre_data.rHx1;
  rHy1            = tyre_data.rHy1;
  rVy1            = tyre_data.rVy1;
  rVy2            = tyre_data.rVy2;
  rVy3            = tyre_data.rVy3;
  rVy4            = tyre_data.rVy4;
  rVy5            = tyre_data.rVy5;
  rVy6            = tyre_data.rVy6;
  lambda__Fz0     = tyre_data.lambda__Fz0;
  lambda__Vyk     = tyre_data.lambda__Vyk;
  lambda__xa      = tyre_data.lambda__xa;
  lambda__yk      = tyre_data.lambda__yk;
  
  [alpha__y, By, Cy, Dy, Ey, SVy, Kya, SHy, mu__y] = MF96_FY0_coeffs(kappa, alpha, phi, Fz, tyre_data);
  

 % main code

  Fz01 = (lambda__Fz0 * Fz0);
  dfz = Fz / Fz01 - 1;
  SHxa = rHx1;
  Bxa = rBx1 * (kappa ^ 2 * rBx2 ^ 2 + 1) ^ (-0.1e1 / 0.2e1) * lambda__xa;
  Cxa = rCx1;
  Dxa = 0.1e1 / cos(Cxa * atan((Bxa * SHxa)));
  Gxa = Dxa * cos(Cxa * atan((Bxa * (alpha + SHxa))));
  SHyk = rHy1;
  DVyk = mu__y * Fz * (dfz * rVy2 + phi * rVy3 + rVy1) * (alpha ^ 2 * rVy4 ^ 2 + 1) ^ (-0.1e1 / 0.2e1);
  SVyk = DVyk * sin(rVy5 * atan((rVy6 * kappa))) * lambda__Vyk;
  Byk = rBy1 * (1 + rBy2 ^ 2 * (alpha - rBy3) ^ 2) ^ (-0.1e1 / 0.2e1) * lambda__yk;
  Cyk = rCy1;
  Dyk = 0.1e1 / cos(Cyk * atan((Byk * SHyk)));
  Gyk = Dyk * cos(Cyk * atan((Byk * (kappa + SHyk))));
  
 end
