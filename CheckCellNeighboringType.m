function neighboring_type = CheckCellNeighboringType(kx_cell_index,
						     ky_cell_index, mesh_cells)
  ## CheckCellNeighboringType - Check the neighboring type for two cells.
  ## @param kx_cell_index The cell index for \f$K_x\f$.
  ## @param ky_cell_index The cell index for \f$K_y\f$.
  ## @param mesh_cells All cells in the mesh.
  ## @param neighboring_type The returned neighboring type for the two cells:
  ## 1: same panel
  ## 2: common edge
  ## 3: common vertex
  ## 4: regular near
  ## 5: regular far

  if (kx_cell_index == ky_cell_index)
    neighboring_type = 1;
  else
    ## Extract the cell node indices from the mesh.
    kx_cell_node_coord_list = mesh_cells(kx_cell_index, :);
    ky_cell_node_coord_list = mesh_cells(ky_cell_index, :);
  endif
endfunction
