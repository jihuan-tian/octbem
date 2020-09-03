function ret = LaplaceHyperSingKernel2DFlat(x, y, nx, ny)
  ## LaplaceHyperSingKernel2DFlat - Hyper singular kernel of Laplace operator in 2D.
  ## @param x A list of source points with dimension N*2.
  ## @param y A list of field points with dimension N*2.
  ## @param nx Normal vector of the cell.
  ## @param ny Normal vector of the cell.

  number_of_points = size(x, 1);
  if (number_of_points != size(y, 1))
    error("The total numbers of the source and field points as well as corresponding normal vectors should be the same!");
  else    
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      dist = norm(x(m, :) - y(m, :));
      ret(m) = 1.0 / 2.0 / pi * (inner_prod(NormalizeVect(nx), NormalizeVect(ny)) / dist^2 - 2 * inner_prod(NormalizeVect(nx), y(m, :) - x(m, :)) * inner_prod(NormalizeVect(ny), y(m, :) - x(m, :)) / dist^4);
    endfor
  endif
endfunction
