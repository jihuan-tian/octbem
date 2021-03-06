function surface_normals = SurfaceNormalOnS2Tria(area_coordinates, node_global_coords, shape_functions_for_geometry, sphere_center, sphere_radius)
  ## SurfaceNormalOnS2Tria - Calculate the surface normal vectors located at
  ## each specified area coordinates on a triangle in S2 manifold.
  ## @param  A list of area coordinates on each of which the global Jacobi
  ## matrix is to be evaluated. They are stored as an N*3 or N*2 matrix
  ## depending on if 3-component or 2-component area coordinates are used.
  ## @param node_global_coords A list of global coordinates for the supporting
  ## nodes in the reference cell in the same order as the shape functions, which
  ## is stored as an N*3 matrix.
  ## @param shape_functions_for_geometry An array of shape function handles used for
  ## describing cell geometry, which are associated with supporting nodes in the
  ## cell.
  ## @param sphere_center Center coordinates of the sphere manifold,
  ## which is stored in a row vector.
  ## @param sphere_radius The radius of the sphere model.
  ## @param surface_normals Normalized surface normal vectors.

  ## Number of nodes in the triangle cell.
  M = size(node_global_coords, 1);
  ## The geometry representation order of the triangle.
  basis_function_order_for_triangle = TriaGeometryOrder(M);

  ## Number of area coordinates and components.
  [number_of_points, N] = size(area_coordinates);
  ## Calculate the functor for calculating the Jacobi matrix of shape functions
  ## with respect to area coordinates. The number of adopted area coordinate
  ## components should be the same as that of the provided area_coordinates.
  switch(N)
    case 2
      shape_function_jacobian_matrix = ShapeFunctionJacobiOn3DTria2Args(basis_function_order_for_triangle);
    case 3
      shape_function_jacobian_matrix = ShapeFunctionJacobiOn3DTria3Args(basis_function_order_for_triangle);
  endswitch

  ## Calculate the surface normal vectors, whose coordinate components are
  ## component Jacobi determinants that are cyclically related to the 3D global
  ## coordinate components, namely, yz, zx and xy.
  surface_normals = zeros(number_of_points, 3);
  
  for n = 1:number_of_points
    jacobian_matrix = PlanarToS2TriaJ(area_coordinates(n, :), shape_functions_for_geometry, node_global_coords, sphere_center, sphere_radius){1} * GlobalJacobiMatrix(area_coordinates(n, :), node_global_coords, shape_function_jacobian_matrix){1};
    
    surface_normals(n, 1) = det(jacobian_matrix([2,3],:));
    surface_normals(n, 2) = det(jacobian_matrix([3,1],:));
    surface_normals(n, 3) = det(jacobian_matrix([1,2],:));
  endfor
  
  ## Normalize the normal vectors.
  surface_normals = surface_normals ./ sqrt(sumsq(surface_normals, 2));
endfunction
