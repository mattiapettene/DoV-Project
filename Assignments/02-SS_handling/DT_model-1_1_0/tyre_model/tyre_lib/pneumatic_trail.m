% Implement the Pneumatic Trail formula
function [trail] = pneumatic_trail(alpha__t, alpha, B, C, D, E)

 % precode

  

 % main code

  t1 = B * alpha__t;
  t2 = atan(t1);
  t6 = atan(-t1 + (t1 - t2) * E);
  t8 = cos(t6 * C);
  t10 = cos(alpha);
  trail = t10 * t8 * D;
  
 end
