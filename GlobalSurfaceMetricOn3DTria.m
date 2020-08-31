function surface_metric = GlobalSurfaceMetricOn3DTria(area_coordinates, node_global_coords)
  ## GlobalSurfaceMetricOn3DTria - Calculate the global surface metric from area
  ## coordinate to global coordinate on a triangle in 3D space.
  ## @param area_coordinates A list of area coordinates on each of
  ## which the global Jacobian.
  ## matrix is to be evaluated. They are stored as an N*3 or N*2 matrix
  ## depending on if 3-component or 2-component area coordinates are used.
  ## @param node_global_coords A list of global coordinates for the supporting
  ## nodes in the reference cell in the same order as the shape functions, which
  ## is stored as an N*3 matrix.
  ## @param surface_metric The calculated surface metric.

  ## Number of nodes in the triangle cell.
  M = size(node_global_coords, 1);
  ## The geometry representation order of the triangle.
  basis_function_order_for_triangle = TriaGeometryOrder(M);

  ## Number of area coordinate components.
  N = size(area_coordinates, 2);
  ## Calculate the functor for calculating the Jacobi matrix of shape functions
  ## with respect to area coordinates. The number of adopted area coordinate
  ## components should be the same as that of the provided area_coordinates.
  switch(N)
    case 2
      shape_function_jacobi_matrix = ShapeFunctionJacobiOn3DTria2Args(basis_function_order_for_triangle);
    case 3
      shape_function_jacobi_matrix = ShapeFunctionJacobiOn3DTria3Args(basis_function_order_for_triangle);
  endswitch

  ## Calculate all the component Jacobian determinants, which are cyclically
  ## related to the 3D global coordinate components, namely, yz, zx and xy.
  jacobi_det_yz = GlobalJacobiDetOn3DTria(area_coordinates, [node_global_coords(:, 2), node_global_coords(:, 3)], shape_function_jacobi_matrix);
  jacobi_det_zx = GlobalJacobiDetOn3DTria(area_coordinates, [node_global_coords(:, 3), node_global_coords(:, 1)], shape_function_jacobi_matrix);
  jacobi_det_xy = GlobalJacobiDetOn3DTria(area_coordinates, [node_global_coords(:, 1), node_global_coords(:, 2)], shape_function_jacobi_matrix);

  ## Calculate the RMS of the componnet Jacobi determinants.
  surface_metric = sqrt(jacobi_det_yz.^2 + jacobi_det_zx.^2 + jacobi_det_xy.^2);
endfunction
