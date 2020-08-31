function h = Tria3DRange(tria_corners)
  ## Tria3DRange - Calculate the range of the triangle in 3D space, which is the
  ## diameter of the circumference circle of the triangle.
  ## @param tria_corners The 3D coordinates for triangle corners, stored as a 3*3 matrix.
  ## @param h The triangle range.

  h = norm(Tria3DCircumcenter(tria_corners) - tria_corners(1, :)) * 2;
endfunction
