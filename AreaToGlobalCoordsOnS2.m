function global_coords = AreaToGlobalCoordsOnS2(area_coord_list, shape_functions_for_geometry, cell_node_coord_list, sphere_center, sphere_radius)
  ## AreaToGlobalCoordsOnS2 - Map a list of area coordinates stored in matrix to
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
  ## @param sphere_center Center coordinates of the sphere manifold,
  ## which is stored in a row vector.
  ## @param sphere_radius The radius of the sphere model.
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
	shape_function_evaluations = zeros(number_of_cell_nodes, 1);
	for n = 1:number_of_cell_nodes
	  shape_function_evaluations(n) = shape_functions_for_geometry{n}(area_coord_list(m, :));
	endfor

	global_coords(m, :) = inner_prod(cell_node_coord_list, shape_function_evaluations);
      endfor

      ## Transform the global coordinates from planar triangle to
      ## curved triangle on S2.
      global_coords = PlanarToS2Tria(global_coords, sphere_center, sphere_radius);
    else
      global_coords = [];
    endif
  endif
endfunction
