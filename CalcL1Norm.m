function L1_norm = CalcL1Norm(dofs, ansatz_function_space, cell_indices_for_integration, mesh_obj, quad_order)
  ## CalcL1Norm - Calculate the L1-norm of the given DOFs on the grid.
  ## @param dofs An array of all DOF values associated to the list of nodes.
  ## @param ansatz_function_space The ansatz function space in a cell,
  ## which contains a list of basis functions.
  ## @param cell_indices_for_integration The list of cell indices for
  ## carrying out the integration.
  ## @param mesh_obj The object holding all the mesh data.
  ## @param quad_order Quadrature order for integration on the cell.
  ## @param L1_norm Calculated L1 norm.

  number_of_cells_for_integration = length(cell_indices_for_integration);
  number_of_dofs_in_cell = length(ansatz_function_space);

  ## Generate quadrature points and weights.
  [qpts_xi1, qpts_xi2, qwts] = triangle_unit_product_set(quad_order);
  qpts = [qpts_xi1', qpts_xi2'];
  qwts = qwts';
  quad_point_num = size(qwts, 1);

  L1_norm = 0;

  ## Iterate over each cell for the integration.
  for m = 1:number_of_cells_for_integration
    ## Current cell index
    e = cell_indices_for_integration(m);
    ## Get dof values in the current cell.
    dofs_in_cell = dofs(mesh_obj.mesh_cells(e, :));
      
    ## Iterate over each quadrature point.
    for q = 1:quad_point_num
      ## The product of the Jacobian and the quadrature weights.
      JxW = qwts(q) * mesh_obj.cell_surface_metrics(e);
      ## Accumulate contribution from each ansatz function.
      local_value = 0;
      for n = 1:number_of_dofs_in_cell
	local_value += dofs_in_cell(n) * ansatz_function_space{n}(qpts(q, :));
      endfor
      L1_norm += abs(local_value) * JxW;
    endfor
  endfor

  L1_norm = L1_norm * triangle_unit_volume();
endfunction
