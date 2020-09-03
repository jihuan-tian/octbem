function ret = LaplaceHyperSingKernel3DFlat(x, y, nx, ny)
  ## LaplaceHyperSingKernel3DFlat - Hyper singular kernel of Laplace operator in 3D.
  ## @param x A list of source points with dimension N*3.
  ## @param y A list of field points with dimension N*3.
  ## @param nx Normal vector of the cell.
  ## @param ny Normal vector of the cell.
  ## @param ret A list of kernel function values at each area coordinate.

  number_of_points = size(x, 1);
  if (number_of_points != size(y, 1))
    error("The total numbers of the source and field points as well as corresponding normal vectors should be the same!");
  else    
    ret = zeros(number_of_points, 1);
    for m = 1:number_of_points
      dist = norm(x(m, :) - y(m, :));
      ret(m) = 1.0 / 4.0 / pi * (inner_prod(NormalizeVect(nx), NormalizeVect(ny)) / dist^3 - 3 * inner_prod(NormalizeVect(nx), x(m, :) - y(m, :)) * inner_prod(NormalizeVect(ny), x(m, :) - y(m, :)) / dist^5);
    endfor
  endif
endfunction
