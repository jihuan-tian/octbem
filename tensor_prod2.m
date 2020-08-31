function ret = tensor_prod2(x, y)
  ## tensor_prod2 - Perform tensor product of two vectors. Note: the second
  ## component varies faster than the first.
  ## @param x The first vector.
  ## @param y The second vector.
  ## @param ret The returned matrix.

  M = length(x);
  N = length(y);
  ret = zeros(M * N, 2);

  counter = 1;
  for m = 1:M
    for n = 1:N
      ret(counter, 1) = x(m);
      ret(counter, 2) = y(n);
      
      counter = counter + 1;
    endfor
  endfor
endfunction
