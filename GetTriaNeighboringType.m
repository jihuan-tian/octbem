function cell_neighboring_type = GetTriaNeighboringType(kx_cell_index, ky_cell_index, mesh_cells)
  ## GetTriaNeighboringType - Get the neighboring type between two cells,
  ## which can be identical, common edge, common vertex or regular.
  ## @param kx_cell_index The cell index for \f$K_x\f$.
  ## @param ky_cell_index The cell index for \f$K_y\f$.
  ## @param mesh_cells All cells in the mesh, each of which stores a list of
  ## ordered node indices.
  ## @param cell_neighboring_type The neighboring type for the two given cells,
  ## which can be:
  ## 1: identical
  ## 2: common edge
  ## 3: common vertex
  ## 4: regular

  if (kx_cell_index == ky_cell_index)
    cell_neighboring_type = 1;
  else
    ## Extract the node indices contained in the two cells.
    kx_cell_node_indices = mesh_cells(kx_cell_index, :);
    ky_cell_node_indices = mesh_cells(ky_cell_index, :);
    ## Determine the geometric order of the two cells.
    kx_cell_order = TriaGeometryOrder(length(kx_cell_node_indices));
    ky_cell_order = TriaGeometryOrder(length(ky_cell_node_indices));
    ## Get corner node indices for the two cells.
    kx_corner_node_indices = kx_cell_node_indices(GetTriaCornerNodeIndices(kx_cell_order));
    ky_corner_node_indices = ky_cell_node_indices(GetTriaCornerNodeIndices(ky_cell_order));

    ## Calculate the intersection of the two node index sets.
    node_index_intersection = intersect(kx_corner_node_indices, ky_corner_node_indices);
    switch(length(node_index_intersection))
      case 0			# Regular
	cell_neighboring_type = 4;
      case 1			# Common vertex
	cell_neighboring_type = 3;
      case 2			# Common edge
	cell_neighboring_type = 2;
      otherwise
	error(cstrcat("Neighbor type detection for the two cells ", num2str(kx_cell_index)," and ", num2str(ky_cell_index), " has an error!"));
    endswitch
  endif
endfunction
