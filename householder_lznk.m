function X = multiplyByQ(X, B, G)
  # B - matrix where vectors are sequential Householder vectors
  # G - matrix where vectors are sequential gamma parameters in H_n
  # returned X = H_1 * .. * H_n * X
  for j = fliplr(1:size(B, 2))
    v = B(:, j);
    g = G(j);
    
    X = X - v * (v' * X) ./ g;
  endfor
endfunction

function [x, R, B, G] = Householder(A, y)
  # variables analogical to http://wazniak.mimuw.edu.pl/index.php?title=MN12
  
  # compute Householder vectors (add 2th norm to nst element of nth vector)
  # and fill ith where i < n with 0
  B = A;
  G = [];
  R = A;
  c = y;
  for j = 1:size(B, 2)
    # create h_n vector (+ if v(1) positive, - otherwise)
    B(:, j) = R(:, j);
    B(1:j-1, j) = 0;
    v1 = B(j, j);
    norm2 = norm(B(:, j), 2);
    signv1 = sign(v1);
    if signv1 == 0
      # case when v1 == 0
      signv1 = 1;
    endif
    B(j, j) = v1 + signv1 * norm2;
    v = B(:, j);
    g = norm2^2 + signv1 * v1 * norm2;
    
    R = R - v * (v' * R) / g;
    c = c - v * (v' * c) / g;
    
    G = [G; g];
  endfor
  # minimalized c = Q^-1 * y and R
  R = R(1:size(A, 2), :);
  c = c(1:size(A, 2), :);
  x = R^-1 * c;
endfunction