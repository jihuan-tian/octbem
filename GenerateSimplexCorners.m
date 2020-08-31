function simplex_corners = GenerateSimplexCorners(N)
  ## GenerateSimplexCorners - Generate a list of corners nodes for a unit
  ## simplex in N-dimensional space.
  ## @param N The dimension of the simplex.
  ## @param simplex_corners The generated simplex corners with a dimension (N+1)*N.

  simplex_corners = zeros(N + 1, N);
  for m = 2:(N + 1)
    simplex_corners(m, m - 1) = 1;
  endfor
endfunction
