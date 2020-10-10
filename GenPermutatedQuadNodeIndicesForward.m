function permutated_node_indices = GenPermutatedQuadNodeIndicesForward(corner_node_no_to_start, std_quad_node_indices)
  ## GenPermutatedQuadNodeIndicesForward - Generate the permutated indices for
  ## the list of nodes in a quadrangle cell in forward zigzag direction by
  ## specifying from which corner node the quadrangle reference cell is
  ## constructed.
  ## @param corner_node_no_to_start The No. for the corner node from which the
  ## reference cell is constructed, which can only be 1, 2, 3 or 4.
  ## @param std_quad_node_indices The standard sequence of node indices in the
  ## reference quadrangle which is stored in a matrix form.
  ## @param permutated_node_indices The returned permutated indices
  ## for the nodes in the quadrangular cell.

  node_index_matrix_dim = size(std_quad_node_indices, 1);
  quad_order = node_index_matrix_dim - 1;
  permutated_node_indices = zeros(1, NumberOfNodesOnQuadrangle(quad_order));

  node_index_counter = 1;
  switch(corner_node_no_to_start)
    case 1
      for n = 1:node_index_matrix_dim
      	for m = node_index_matrix_dim:(-1):1
      	  permutated_node_indices(node_index_counter) = std_quad_node_indices(m, n);
      	  node_index_counter = node_index_counter + 1;
      	endfor
      endfor
    case 2
      for m = 1:node_index_matrix_dim
	for n = 1:node_index_matrix_dim
	  permutated_node_indices(node_index_counter) = std_quad_node_indices(m, n);
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
    case 3
      for m = node_index_matrix_dim:(-1):1
	for n = node_index_matrix_dim:(-1):1
	  permutated_node_indices(node_index_counter) = std_quad_node_indices(m, n);
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
    case 4
      for n = node_index_matrix_dim:(-1):1
	for m = 1:node_index_matrix_dim
	  permutated_node_indices(node_index_counter) = std_quad_node_indices(m, n);
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
  endswitch
endfunction
