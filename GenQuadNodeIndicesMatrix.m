function quad_node_indices_matrix = GenQuadNodeIndicesMatrix(quad_node_indices)
  ## GenQuadNodeIndicesMatrix - Generate the rectangle for the list of node
  ## indices in a quadrangle.
  ## @param quad_node_indices A list of node indices for the triangle.
  ## @param quad_node_indices_matrix The generated node indices matrix.

  number_of_nodes = length(quad_node_indices);
  quad_geometry_order = QuadGeometryOrder(number_of_nodes);
  
  quad_node_indices_matrix_dim = quad_geometry_order + 1;
  quad_node_indices_matrix = zeros(quad_node_indices_matrix_dim, quad_node_indices_matrix_dim);

  node_index_counter = 1;
  for n = 1:quad_node_indices_matrix_dim
    for m = quad_node_indices_matrix_dim:(-1):1
      quad_node_indices_matrix(m, n) = quad_node_indices(node_index_counter);
      node_index_counter = node_index_counter + 1;
    endfor
  endfor
endfunction
