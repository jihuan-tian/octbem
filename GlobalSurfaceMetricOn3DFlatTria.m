function surface_metric = GlobalSurfaceMetricOn3DFlatTria(corner_global_coords)
  ## GlobalSurfaceMetricOn3DFlatTria - Calculate the global surface metric on a
  ## flat triangle in 3D space, which is the double of triangle area.
  ## @param corner_global_coords A list of global coordinates for the corner
  ## nodes in the reference cell.
  ## @param surface_metric The calculated surface metric, which is
  ## ensured to be positive.

  edge1 = corner_global_coords(2, :) - corner_global_coords(1, :);
  edge2 = corner_global_coords(3, :) - corner_global_coords(2, :);

  surface_metric = norm(cross(edge1, edge2));
endfunction
