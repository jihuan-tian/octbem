function ret = tensor_prod4(w, x, y, z)
  ## tensor_prod4 - Perform tensor product of four vectors. Note: the fourth
  ## component runs the fastest, then the third, then the second and then the
  ## first.
  ## @param w The first vector.
  ## @param x The second vector.
  ## @param y The third vector.
  ## @param z The fourth vector.
  ## @param ret The returned matrix.

  K = length(w);
  L = length(x);
  M = length(y);
  N = length(z);
  ret = zeros(K * L * M * N, 4);

  counter = 1;
  for k = 1:K
    for l = 1:L
      for m = 1:M
	for n = 1:N
	  ret(counter, 1) = w(k);
	  ret(counter, 2) = x(l);
	  ret(counter, 3) = y(m);
	  ret(counter, 4) = z(n);
	  
	  counter = counter + 1;
	endfor
      endfor
    endfor
  endfor
endfunction
