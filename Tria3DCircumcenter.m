function circumcenter = Tria3DCircumcenter(tria_corners)
  ## Tria3DCircumcenter - Calculate the circumcenter for a triangle in 3D space.
  ## The calculation formula is:
  ##           |c-a|^2 [(b-a)x(c-a)]x(b-a) + |b-a|^2 (c-a)x[(b-a)x(c-a)]
  ## m = a + ---------------------------------------------------------.
  ##                            2 | (b-a)x(c-a) |^2
  ## @param tria_corners The 3D coordinates for triangle corners, stored as a 3*3 matrix.
  ## @param circumcenter_global_coord The calculated 3D coordinate for the circumcenter.

  ab = tria_corners(2, :) - tria_corners(1, :);
  ac = tria_corners(3, :) - tria_corners(1, :);

  circumcenter = tria_corners(1, :) + (norm(ac)^2 * cross(cross(ab, ac), ab) + norm(ab)^2 * cross(ac, cross(ab, ac))) / (2 * norm(cross(ab, ac))^2);
endfunction
