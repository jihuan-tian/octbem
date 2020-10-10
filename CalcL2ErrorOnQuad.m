function L2_error = CalcL2ErrorOnQuad(dofs, functor, ansatz_function_space, shape_functions_for_geometry, cell_indices_for_integration, mesh_obj, qpts, qwts)
  ## CalcL2ErrorOnQuad - Calculate the L1-error between a grid
  ## function and an analytical expression for quadrangle mesh.
  ## @param dofs An array of all DOF values associated to the list of nodes.
  ## @param functor The functor defining the analytical expression
  ## dependent on the global coordinates.
  ## @param ansatz_function_space The ansatz function space in a cell,
  ## which contains a list of basis functions.
  ## @param shape_functions_for_geometry The shape function space for
  ## describing the cell geometry.
  ## @param cell_indices_for_integration The list of cell indices for
  ## carrying out the integration.
  ## @param mesh_obj The object holding all the mesh data.
  ## @param qpts Quadrature points on the quadrangle.
  ## @param qwts Quadrature weights on the quadrature.
  ## @param L2_error Calculated relative L2 error.

  number_of_cells_for_integration = length(cell_indices_for_integration);
  number_of_dofs_in_cell = length(ansatz_function_space);

  quad_point_num = size(qwts, 1);

  ## The L2 norm of the difference between the grid function and the
  ## analytical expression.
  diff_L2_norm = 0;

  ## Iterate over each cell for the integration.
  for m = 1:number_of_cells_for_integration
    ## Current cell index
    e = cell_indices_for_integration(m);
    ## Extract cell node coordinates.
    cell_node_coord_list = mesh_obj.mesh_nodes(mesh_obj.mesh_cells(e, :), :);
    ## Get dof values in the current cell.
    dofs_in_cell = dofs(mesh_obj.mesh_cells(e, :));

    ## Construct surface metric functor on the current cell.
    J = @(area_coord) GlobalSurfaceMetricOn3DQuad(area_coord, cell_node_coord_list);
      
    ## Iterate over each quadrature point.
    for q = 1:quad_point_num
      ## The product of the Jacobian and the quadrature weights.
      JxW = qwts(q) * J(qpts(q, :));
      ## Accumulate contribution from each ansatz function for
      ## evaluating the grid function at the quadrature point.
      local_value = 0;
      for n = 1:number_of_dofs_in_cell
	local_value += dofs_in_cell(n) * ansatz_function_space{n}(qpts(q, :));
      endfor
      ## Get the global coordinates for the current quadrature point
      ## and evaluate the functor at this point.
      global_coords = AreaToGlobalCoords(qpts(q, :), shape_functions_for_geometry, cell_node_coord_list);
      diff_L2_norm += (local_value - functor(global_coords))^2 * JxW;
    endfor
  endfor

  diff_L2_norm = sqrt(diff_L2_norm);

  ## Calculate the L2 norm of the functor.
  functor_L2_norm = CalcL2NormForExprOnQuad(functor, shape_functions_for_geometry, cell_indices_for_integration, mesh_obj, qpts, qwts);

  L2_error = diff_L2_norm / functor_L2_norm;
endfunction
