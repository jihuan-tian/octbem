function tria_node_indices_matrix = GenTriaNodeIndicesMatrix(tria_node_indices)
  ## GenTriaNodeIndicesMatrix - Generate the lower triangle for the list of node
  ## indices in a triangle.
  ## @param tria_node_indices A list of node indices for the triangle.
  ## @param tria_node_indices_matrix The generated node indices matrix.

  number_of_nodes = length(tria_node_indices);
  tria_geometry_order = TriaGeometryOrder(number_of_nodes);
  
  tria_node_indices_matrix_dim = tria_geometry_order + 1;
  tria_node_indices_matrix = zeros(tria_node_indices_matrix_dim, tria_node_indices_matrix_dim);

  node_index_counter = 1;
  for n = 1:tria_node_indices_matrix_dim
    for m = tria_node_indices_matrix_dim:(-1):n
      tria_node_indices_matrix(m, n) = tria_node_indices(node_index_counter);
      node_index_counter = node_index_counter + 1;
    endfor
  endfor
endfunction
