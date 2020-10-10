function surface_normal = SurfaceNormalOn3DFlatTria(corner_global_coords)
  ## SurfaceNormalOn3DFlatTria - Calculate the surface normal vectors located at
  ## each specified area coordinates for a flat triangle.
  ## @param corner_global_coords A list of global coordinates for the corner
  ## nodes in the reference cell.
  ## @param surface_normal The calculated surface normal vector. Its
  ## direction is the opposite of the original mesh direction.

  edge1 = corner_global_coords(2, :) - corner_global_coords(1, :);
  edge2 = corner_global_coords(3, :) - corner_global_coords(2, :);
  
  surface_normal = NormalizeVect(cross(edge2, edge1));
endfunction
