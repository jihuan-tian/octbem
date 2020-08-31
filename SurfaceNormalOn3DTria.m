function surface_normals = SurfaceNormalOn3DTria(area_coordinates, node_global_coords)
  ## SurfaceNormalOn3DTria - Calculate the surface normal vectors located at
  ## each specified area coordinates.
  ## @param  A list of area coordinates on each of which the global Jacobi
  ## matrix is to be evaluated. They are stored as an N*3 or N*2 matrix
  ## depending on if 3-component or 2-component area coordinates are used.
  ## @param node_global_coords A list of global coordinates for the supporting
  ## nodes in the reference cell in the same order as the shape functions, which
  ## is stored as an N*3 matrix.
  ## @param surface_normals The calculated surface normal vectors (normalized).

  ## Number of nodes in the triangle cell.
  M = size(node_global_coords, 1);
  ## The geometry representation order of the triangle.
  basis_function_order_for_triangle = TriaGeometryOrder(M);

  ## Number of area coordinates and components.
  [L, N] = size(area_coordinates);
  ## Calculate the functor for calculating the Jacobi matrix of shape functions
  ## with respect to area coordinates. The number of adopted area coordinate
  ## components should be the same as that of the provided area_coordinates.
  switch(N)
    case 2
      shape_function_jacobi_matrix = ShapeFunctionJacobiOn3DTria2Args(basis_function_order_for_triangle);
    case 3
      shape_function_jacobi_matrix = ShapeFunctionJacobiOn3DTria3Args(basis_function_order_for_triangle);
  endswitch

  ## Calculate the surface normal vectors, whose coordinate components are
  ## component Jacobi determinants that are cyclically related to the 3D global
  ## coordinate components, namely, yz, zx and xy.
  surface_normals = zeros(L, 3);
  surface_normals(:, 1) = GlobalJacobiDetOn3DTria(area_coordinates, [node_global_coords(:, 2), node_global_coords(:, 3)], shape_function_jacobi_matrix);
  surface_normals(:, 2) = GlobalJacobiDetOn3DTria(area_coordinates, [node_global_coords(:, 3), node_global_coords(:, 1)], shape_function_jacobi_matrix);
  surface_normals(:, 3) = GlobalJacobiDetOn3DTria(area_coordinates, [node_global_coords(:, 1), node_global_coords(:, 2)], shape_function_jacobi_matrix);
  ## Normalize the normal vectors.
  surface_normals = surface_normals ./ sqrt(sumsq(surface_normals, 2));
endfunction
