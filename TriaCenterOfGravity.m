function cog = TriaCenterOfGravity(tria_corners)
  ## TriaCenterOfGravity - Calculte the coordinates of center of
  ## gravity for the triangle.
  ## @param tria_corners The coordinates for triangle corners, stored as a 3*d matrix, where d is the space dimension.
  ## @param cog The calculated center of gravity.

  space_dim = size(tria_corners, 2);
  cog = zeros(1, space_dim);
  
  for d = 1:space_dim
    cog(d) = average(tria_corners(:, d));
  endfor
endfunction
