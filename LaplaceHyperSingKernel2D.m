function ret = LaplaceHyperSingKernel2D(x, y, nx, ny)
  ## LaplaceHyperSingKernel2D - Hyper singular kernel of Laplace operator in 2D.
  ## @param x A list of source points with dimension N*2.
  ## @param y A list of field points with dimension N*2.
  ## @param nx A list of normal vectors at each \f$x\f$.
  ## @param ny A list of normal vectors at each \f$y\f$.

  number_of_points = size(x, 1);
  if ((number_of_points != size(y, 1)) || 
      (number_of_points != size(nx, 1)) ||
      (number_of_points != size(ny, 1)))
    error("The total numbers of the source and field points as well as corresponding normal vectors should be the same!");
  else    
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      dist = norm(x(m, :) - y(m, :));
      ret(m) = 1.0 / 2.0 / pi * (inner_prod(NormalizeVect(nx(m, :)), NormalizeVect(ny(m, :))) / dist^2 - 2 * inner_prod(NormalizeVect(nx(m, :)), y(m, :) - x(m, :)) * inner_prod(NormalizeVect(ny(m, :)), y(m, :) - x(m, :)) / dist^4);
    endfor
  endif
endfunction
