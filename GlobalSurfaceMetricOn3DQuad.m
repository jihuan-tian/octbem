function surface_metrics = GlobalSurfaceMetricOn3DQuad(area_coordinates, node_global_coords)
  ## GlobalSurfaceMetricOn3DQuad - Calculate the global surface metric from area
  ## coordinate to global coordinate on a quadrangle in 3D space.
  ## @param area_coordinates A list of area coordinates on each of
  ## which the global Jacobian matrix is to be evaluated. They are
  ## stored as an N*2 matrix.
  ## @param node_global_coords A list of global coordinates for the supporting
  ## nodes in the reference cell in the same order as the shape functions, which
  ## is stored as an N*3 matrix.
  ## @param surface_metrics The calculated surface metric at each point.

  ## Number of nodes in the quadrangle cell.
  M = size(node_global_coords, 1);
  ## The geometry representation order of the quadrangle.
  basis_function_order_for_quadrangle = QuadGeometryOrder(M);


  number_of_points = size(area_coordinates, 1);
  shape_function_jacobi_matrix = ShapeFunctionJacobiOn3DQuad(basis_function_order_for_quadrangle);
  
  ## Calculate all the component Jacobian determinants, which are cyclically
  ## related to the 3D global coordinate components, namely, yz, zx and xy.
  surface_metrics = zeros(number_of_points, 1);

  for n = 1:number_of_points
    jacobi_det_yz = GlobalJacobiDet(area_coordinates(n, :), [node_global_coords(:, 2), node_global_coords(:, 3)], shape_function_jacobi_matrix);
    jacobi_det_zx = GlobalJacobiDet(area_coordinates(n, :), [node_global_coords(:, 3), node_global_coords(:, 1)], shape_function_jacobi_matrix);
    jacobi_det_xy = GlobalJacobiDet(area_coordinates(n, :), [node_global_coords(:, 1), node_global_coords(:, 2)], shape_function_jacobi_matrix);

    ## Calculate the RMS of the componnet Jacobi determinants.
    surface_metrics(n) = sqrt(jacobi_det_yz.^2 + jacobi_det_zx.^2 + jacobi_det_xy.^2);
  endfor
endfunction
