function cog = PolygonCenterOfGravity(polygon_corners)
  ## PolygonCenterOfGravity - Calculte the coordinates of center of
  ## gravity for the polygon.
  ## @param polygon_corners The coordinates for polygon corners,
  ## stored as a n*d matrix, where n is the number of corner points
  ## and d is the space dimension.
  ## @param cog The calculated center of gravity.

  space_dim = size(polygon_corners, 2);
  cog = zeros(1, space_dim);
  
  for d = 1:space_dim
    cog(d) = average(polygon_corners(:, d));
  endfor
endfunction
