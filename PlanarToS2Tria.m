function s2_tria_coords = PlanarToS2Tria(planar_tria_coords,
					 sphere_center, sphere_radius)
  ## PlanarToS2Tria - Transform point coordinates on planar triangle
  ## to triangle on sphere.
  ## @param planar_tria_coords A list of point coordinates on planar
  ## triangle.
  ## @param sphere_center Center coordinates of the sphere manifold,
  ## which is stored in a row vector.
  ## @param sphere_radius The radius of the sphere model.
  ## @param s2_tria_coords A list of transformed point coordinates on
  ## S2 triangle.

  for m = 1:size(planar_tria_coords, 1)
    s2_tria_coords(m, :) = (planar_tria_coords(m, :) - sphere_center) / norm(planar_tria_coords(m, :) - sphere_center) * sphere_radius + sphere_center;
  endfor
endfunction
