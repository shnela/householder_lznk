WITHOUT_IT_WNOT_WORK = 1; # DAFAQ

function [x, R, B] = Householder(A, y)
  # compute Householder vectors
  B = A;
  B
  for v = B
    #TODO ifelse +/-
    rank(v, 2)
    v(1) += rank(v, 2);
  endfor
  B
  
  x = A \ y; # TODO
  [_, R] = qr(A); #TODO
endfunction