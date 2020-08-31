function ret = LaplaceSLPKernel3D(x, y)
  ## LaplaceSLPKernel3D - Integration kernel for the single layer potential of
  ## Laplace operator in 3D.
  ## @param x A list of source points with dimension N*3.
  ## @param y A list of field points with dimension N*3.

  number_of_points = size(x, 1);
  if (number_of_points != size(y, 1))
    error("The total number of the source points should be the same as field points!");
  else
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      ret(m) = 1.0 / 4.0 / pi / norm(x(m, :) - y(m, :));
    endfor
  endif
endfunction
