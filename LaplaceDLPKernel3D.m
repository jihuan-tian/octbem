function ret = LaplaceDLPKernel3D(x, y, nx, ny)
  ## LaplaceDLPKernel3D - Integration kernel for the double layer potential of
  ## Laplace operator in 3D.
  ## @param x A list of source points with dimension N*3.
  ## @param y A list of field points with dimension N*3.
  ## @param nx A list of normal vectors at each \f$x\f$, but this is just a placeholder and not actually used.
  ## @param ny A list of normal vectors at each \f$y\f$.

  number_of_points = size(x, 1);
  if ((number_of_points != size(y, 1)) || (number_of_points != size(nx, 1)) || (number_of_points != size(ny, 1)))
    error("The total numbers of the source and field points as well as corresponding normal vectors should be the same!");
  else    
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      ret(m) = inner_prod(x(m, :) - y(m, :), NormalizeVect(ny(m, :))) / 4.0 / pi / norm(x(m, :) - y(m, :))^3;
    endfor
  endif
endfunction
