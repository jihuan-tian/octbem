function ret = LaplaceDLPKernel3DFlat(x, y, nx, ny)
  ## LaplaceDLPKernel3DFlat - Integration kernel for the double layer potential of
  ## Laplace operator in 3D.
  ## @param x A list of source points with dimension N*3.
  ## @param y A list of field points with dimension N*3.
  ## @param nx Normal vector of the cell, but this is just a placeholder and not actually used.
  ## @param ny Normal vector of the cell.

  number_of_points = size(x, 1);
  if (number_of_points != size(y, 1))
    error("The total numbers of the source and field points as well as corresponding normal vectors should be the same!");
  else    
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      ret(m) = inner_prod(x(m, :) - y(m, :), NormalizeVect(ny)) / 4.0 / pi / norm(x(m, :) - y(m, :))^3;
    endfor
  endif
endfunction
