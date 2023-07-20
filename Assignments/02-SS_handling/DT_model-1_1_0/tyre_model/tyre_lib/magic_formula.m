% Implement the Magic formula
function [mf_res] = magic_formula(x, B, C, D, E, SV)

 % precode

  

 % main code

  t1 = B * x;
  t2 = atan(t1);
  t6 = atan(-t1 + (t1 - t2) * E);
  t8 = sin(t6 * C);
  mf_res = -t8 * D + SV;
  
 end
