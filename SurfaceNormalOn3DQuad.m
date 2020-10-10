function surface_normals = SurfaceNormalOn3DQuad(area_coordinates, node_global_coords)
  ## SurfaceNormalOn3DQuad - Calculate the surface normal vectors located at
  ## each specified area coordinates.
  ## @param area_coordinates A list of area coordinates on each of which the global Jacobi
  ## matrix is to be evaluated. They are stored as an N*2 matrix
  ## @param node_global_coords A list of global coordinates for the supporting
  ## nodes in the reference cell in the same order as the shape functions, which
  ## is stored as an N*3 matrix.
  ## @param surface_normals The calculated surface normal vectors (normalized).

  ## Number of nodes in the quadrangle cell.
  M = size(node_global_coords, 1);
  ## The geometry representation order of the quadrangle.
  basis_function_order_for_quadrangle = QuadGeometryOrder(M);

  ## Number of area coordinates and components.
  number_of_points = size(area_coordinates, 1);
  shape_function_jacobi_matrix = ShapeFunctionJacobiOn3DQuad(basis_function_order_for_quadrangle);

  ## Calculate the surface normal vectors, whose coordinate components are
  ## component Jacobi determinants that are cyclically related to the 3D global
  ## coordinate components, namely, yz, zx and xy.
  surface_normals = zeros(number_of_points, 3);
  surface_normals(:, 1) = GlobalJacobiDet(area_coordinates, [node_global_coords(:, 2), node_global_coords(:, 3)], shape_function_jacobi_matrix);
  surface_normals(:, 2) = GlobalJacobiDet(area_coordinates, [node_global_coords(:, 3), node_global_coords(:, 1)], shape_function_jacobi_matrix);
  surface_normals(:, 3) = GlobalJacobiDet(area_coordinates, [node_global_coords(:, 1), node_global_coords(:, 2)], shape_function_jacobi_matrix);
  ## Normalize the normal vectors.
  surface_normals = surface_normals ./ sqrt(sumsq(surface_normals, 2));
endfunction
