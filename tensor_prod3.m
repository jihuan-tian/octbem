function ret = tensor_prod3(x, y, z)
  ## tensor_prod3 - Perform tensor product of three vectors. Note: the third
  ## component runs the fastest, then the second, and then the first.
  ## @param x The first vector.
  ## @param y The second vector.
  ## @param z The third vector.
  ## @param ret The returned matrix.

  L = length(x);
  M = length(y);
  N = length(z);
  ret = zeros(L * M * N, 3);

  counter = 1;
  for l = 1:L
    for m = 1:M
      for n = 1:N
	ret(counter, 1) = x(l);
	ret(counter, 2) = y(m);
	ret(counter, 3) = z(n);
	
	counter = counter + 1;
      endfor
    endfor
  endfor
endfunction
