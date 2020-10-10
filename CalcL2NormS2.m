function L2_norm = CalcL2NormS2(dofs, ansatz_function_space, shape_functions_for_geometry, cell_indices_for_integration, mesh_obj, quad_order, sphere_center, sphere_radius)
  ## CalcL2NormS2 - Calculate the L2-norm of the given DOFs on the grid.
  ## @param dofs An array of all DOF values associated to the list of nodes.
  ## @param ansatz_function_space The ansatz function space in a cell,
  ## which contains a list of basis functions.
  ## @param shape_functions_for_geometry The shape function space for
  ## describing the geometry of the cell.
  ## @param cell_indices_for_integration The list of cell indices for
  ## carrying out the integration.
  ## @param mesh_obj The object holding all the mesh data.
  ## @param quad_order Quadrature order for integration on the cell.
  ## @param sphere_center Center coordinates of the sphere manifold,
  ## which is stored in a row vector.
  ## @param sphere_radius The radius of the sphere model.
  ## @param L2_norm Calculated L2 norm.

  number_of_cells_for_integration = length(cell_indices_for_integration);
  number_of_dofs_in_cell = length(ansatz_function_space);

  ## Generate quadrature points and weights.
  [qpts_xi1, qpts_xi2, qwts] = triangle_unit_product_set(quad_order);
  qpts = [qpts_xi1', qpts_xi2'];
  qwts = qwts';
  quad_point_num = size(qwts, 1);

  L2_norm = 0;

  ## Iterate over each cell for the integration.
  for m = 1:number_of_cells_for_integration
    ## Current cell index
    e = cell_indices_for_integration(m);
    ## Get dof values in the current cell.
    dofs_in_cell = dofs(mesh_obj.mesh_cells(e, :));
    ## Generate the surface Jacobian functor.
    J = @(x) GlobalSurfaceMetricOnS2Tria(x, mesh_obj.mesh_nodes(mesh_obj.mesh_cells(e, :), :), shape_functions_for_geometry, sphere_center, sphere_radius);
      
    ## Iterate over each quadrature point.
    for q = 1:quad_point_num
      ## The product of the Jacobian and the quadrature weights.
      JxW = qwts(q) * J(qpts(q, :));
      ## Accumulate contribution from each ansatz function.
      local_value = 0;
      for n = 1:number_of_dofs_in_cell
	local_value += dofs_in_cell(n) * ansatz_function_space{n}(qpts(q, :));
      endfor
      L2_norm += local_value^2 * JxW;
    endfor
  endfor

  L2_norm = sqrt(L2_norm * triangle_unit_volume());
endfunction
