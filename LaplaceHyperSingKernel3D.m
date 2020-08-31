function ret = LaplaceHyperSingKernel3D(x, y, nx, ny)
  ## LaplaceHyperSingKernel3D - Hyper singular kernel of Laplace operator in 3D.
  ## @param x A list of source points with dimension N*3.
  ## @param y A list of field points with dimension N*3.
  ## @param nx A list of normal vectors at each \f$x\f$.
  ## @param ny A list of normal vectors at each \f$y\f$.
  ## @param ret A list of kernel function values at each area coordinate.

  number_of_points = size(x, 1);
  if ((number_of_points != size(y, 1)) || (number_of_points != size(nx, 1)) || (number_of_points != size(ny, 1)))
    error("The total numbers of the source and field points as well as corresponding normal vectors should be the same!");
  else    
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      dist = norm(x(m, :) - y(m, :));
      ret(m) = 1.0 / 4.0 / pi * (inner_prod(NormalizeVect(nx(m, :)), NormalizeVect(ny(m, :))) / dist^3 - 3 * inner_prod(NormalizeVect(nx(m, :)), x(m, :) - y(m, :)) * inner_prod(NormalizeVect(ny(m, :)), x(m, :) - y(m, :)) / dist^5);
    endfor
  endif
endfunction
