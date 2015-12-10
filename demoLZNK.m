source("householder_lznk.m")

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



function runTest(A, y)
  [x, R, B, G] = Householder(A,y);
  
  disp("Test A (coefficient accuracy)")
  x
  
  disp("Run test B (QR factorization accuracy)")
  [m, n] = size(A);
  R(n+1:m, :) = 0;
  QR = multiplyByQ(R, B, G);
  E = norm(A - QR, 2);
  E
endfunction

for m = [10, 20, 100]
  disp("Test data for m = "), disp(m)
  [A, y] = generateAy(m);
  runTest(A, y)
endfor