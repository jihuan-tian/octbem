function global_coords = AreaToGlobalCoords(area_coord_list, shape_functions_for_geometry, cell_node_coord_list)
  ## AreaToGlobalCoords - Map a list of area coordinates stored in matrix to
  ## global coordinates.
  ## @param area_coord_list A list of area coordinates stored in a matrix with
  ## dimension N*M, where N is the number of area coordiantes, M is the number
  ## of coordinate components.
  ## @param shape_functions_for_geometry An array of shape function handles used for
  ## describing cell geometry, which are associated with supporting nodes in the
  ## cell.
  ## @param cell_node_list A list of global coordinates for cell nodes with
  ## dimension N*M, where N is the number of cell nodes, M is the space
  ## dimension.
  ## @param global_coords The transformed global coordinates.


  [number_of_area_coords, area_coord_dim] = size(area_coord_list);
  [number_of_cell_nodes, space_dim] = size(cell_node_coord_list);
  number_of_shape_functions = length(shape_functions_for_geometry);

  if (number_of_cell_nodes != number_of_shape_functions)
    global_coords = [];
    error("The number of shape functions for describing cell geometry should be the same as the number of cell nodes!");
  else
    if (number_of_area_coords > 0)
      global_coords = zeros(number_of_area_coords, space_dim);

      for m = 1:number_of_area_coords
	## Linear combination of each shape functions by weighting them with
	## cell node coordinates.
	for n = 1:number_of_cell_nodes
	  global_coords(m, :) = global_coords(m, :) + cell_node_coord_list(n, :) * shape_functions_for_geometry{n}(area_coord_list(m, :));
	endfor
      endfor
    else
      global_coords = [];
    endif
  endif
endfunction
