WITHOUT_IT_WNOT_WORK = 1; # DAFAQ

function [A, y] = generateAy(m)
  # construct matrix A (x^2 + x + 1) for x in range [1:10]
  A = [];
  for x = linspace(0, 10, m)
    v = [x^2, x, 1];
    A = [A; v];
  endfor
  
  # generate results with perturbations
  e = rand(m, 1) * 2 * 10^-2 - 10^-2;
  y = A * [1; -5; 2] + e;
endfunction



function runTestA(A, y)
  Householder(A,y)
endfunction

function runTestB(A, y)
  [_, R, B] = Householder(A,y);
  #E = A - B * R;
  #norm(E, 2)
endfunction

#for m = [10, 20, 100]
for m = [10]
  [A, y] = generateAy(11);
  runTestA(A, y)
  runTestB(A, y)
endfor