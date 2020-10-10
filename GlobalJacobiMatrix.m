function global_jacobi_matrices = GlobalJacobiMatrix(area_coordinates, node_global_coord_components, shape_function_jacobi_matrix)
  ## GlobalJacobiMatrix - Calculate the Jacobi matrix from area coordinate to
  ## global coordinate on a triangle in 3D space.
  ## @param  A list of area coordinates on each of which the global Jacobi
  ## matrix is to be evaluated. They are stored as an N*3 or N*2 matrix
  ## depending on if 3-component or 2-component area coordinates are used.
  ## @param node_global_coord_components A list of global coordinates with only
  ## the interested components (for example, yz, zx or xy components) for the
  ## supporting nodes in the reference cell in the same order as the shape
  ## functions, which is stored as an N*2 matrix.
  ## @param shape_function_jacobi_matrix The functor for calculating the Jacobi
  ## matrix of shape functions with respect to area coordinates. The number of
  ## adopted area coordinate components should be the same as that of the
  ## provided area_coordinates.
  ## @param global_jacobi_matrices The calculated list of Jacobi matrices
  ## evaluated at each area coordinate.

  ## Number of area coordinates for the evaluation of Jacobi matrix.
  N = size(area_coordinates,  1);
  global_jacobi_matrices = cell(N, 1);

  for m = 1:N
    global_jacobi_matrices{m} = node_global_coord_components' * shape_function_jacobi_matrix(area_coordinates(m, :));
  endfor
endfunction
