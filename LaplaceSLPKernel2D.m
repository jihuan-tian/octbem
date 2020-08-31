function ret = LaplaceSLPKernel2D(x, y)
  ## LaplaceSLPKernel2D - Integration kernel for the single layer potential of
  ## Laplace operator in 2D.
  ## @param x A list of source points with dimension N*2.
  ## @param y A list of field points with dimension N*2.

  number_of_points = size(x, 1);
  if (number_of_points != size(y, 1))
    error("The total number of the source points should be the same as field points!");
  else
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      ret(m) = -1.0 / 2.0 / pi * log(1.0 / norm(x(m, :) - y(m, :)));
    endfor
  endif
endfunction
