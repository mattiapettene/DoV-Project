% Magic formula - stiffness
function [mf_res] = magic_formula_stiffness(x, B, C, D, E, SV)

 % precode

  

 % main code

  t2 = (B ^ 2);
  t3 = (x ^ 2);
  t11 = B * x;
  t12 = atan(t11);
  t15 = -t11 + (t11 - t12) * E;
  t16 = t15 ^ 2;
  t20 = atan(t15);
  t22 = cos(t20 * C);
  mf_res = -t22 / (t16 + 0.1e1) * (-B + (B - 1 / (t3 * t2 + 1) * B) * E) * D * C;
  
 end
