function permutated_node_indices = GenPermutatedTriaNodeIndicesForward(corner_node_no_to_start, std_tria_node_indices)
  ## GenPermutatedTriaNodeIndicesForward - Generate the permutated indices for
  ## the list of nodes in a triangle cell in forward zigzag direction by
  ## specifying from which corner node the triangle reference cell is
  ## constructed.
  ## @param corner_node_no_to_start The No. for the corner node from which the
  ## reference cell is constructed, which can only be 1, 2 or 3.
  ## @param std_tria_node_indices The standard sequence of node indices in the
  ## reference triangle.
  ## @param permutated_node_indices The returned permutated indices for the nodes in the
  ## triangle cell.

  node_index_matrix_dim = size(std_tria_node_indices, 1);
  tria_order = node_index_matrix_dim - 1;
  permutated_node_indices = zeros(1, NumberOfNodesOnTriangle(tria_order));

  node_index_counter = 1;
  switch(corner_node_no_to_start)
    case 1
      for n = 1:node_index_matrix_dim
      	for m = node_index_matrix_dim:(-1):n
      	  permutated_node_indices(node_index_counter) = std_tria_node_indices(m, n);
      	  node_index_counter = node_index_counter + 1;
      	endfor
      endfor
    case 2
      for n = 1:node_index_matrix_dim
	col_idx = 1;
	for m = n:node_index_matrix_dim
	  permutated_node_indices(node_index_counter) = std_tria_node_indices(m, col_idx);
	  col_idx = col_idx + 1;
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
    case 3
      for m = node_index_matrix_dim:(-1):1
	for n = m:(-1):1
	  permutated_node_indices(node_index_counter) = std_tria_node_indices(m, n);
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
  endswitch
endfunction
