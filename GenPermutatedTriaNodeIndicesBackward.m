function permutated_node_indices = GenPermutatedTriaNodeIndicesBackward(corner_node_no_to_start, tria_node_indices_matrix)
  ## GenPermutatedTriaNodeIndicesBackward - Generate the permutated indices for
  ## the list of nodes in a triangle cell in backward zigzag direction by
  ## specifying from which corner node the triangle reference cell is
  ## constructed.
  ## @param corner_node_no_to_start The No. for the corner node from which the
  ## reference cell is constructed, which can only be 1, 2 or 3.
  ## @param tria_node_indices_matrix The standard sequence of node indices in the
  ## reference triangle.
  ## @param permutated_node_indices The returned permutated indices for the nodes in the
  ## triangle cell.

  node_index_matrix_dim = size(tria_node_indices_matrix, 1);
  tria_order = node_index_matrix_dim - 1;
  permutated_node_indices = zeros(1, NumberOfNodesOnTriangle(tria_order));

  node_index_counter = 1;
  switch(corner_node_no_to_start)
    case 1
      for m = node_index_matrix_dim:(-1):1
	for n = 1:m
	  permutated_node_indices(node_index_counter) = tria_node_indices_matrix(m, n);
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
    case 2
      for n = 1:node_index_matrix_dim
	for m = n:node_index_matrix_dim
	  permutated_node_indices(node_index_counter) = tria_node_indices_matrix(m, n);
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
    case 3
      for m = node_index_matrix_dim:(-1):1
	row_idx = node_index_matrix_dim;
	for n = m:(-1):1
	  permutated_node_indices(node_index_counter) = tria_node_indices_matrix(row_idx, n);
	  row_idx = row_idx - 1;
	  node_index_counter = node_index_counter + 1;
	endfor
      endfor
  endswitch
endfunction
