function L1_norm = CalcL1NormForExpr(functor, shape_functions_for_geometry,
				     cell_indices_for_integration, mesh_obj, quad_order)
  ## CalcL1NormForExpr - Calculate the L1-norm of the given analytical
  ## expression.
  ## @param functor The functor defining the analytical expression
  ## dependent on the global coordinates.
  ## @param shape_functions_for_geometry The shape function space for
  ## describing the cell geometry.
  ## @param cell_indices_for_integration The list of cell indices for
  ## carrying out the integration.
  ## @param mesh_obj The object holding all the mesh data.
  ## @param quad_order Quadrature order for integration on the cell.
  ## @param L1_norm Calculated L1 norm.

  number_of_cells_for_integration = length(cell_indices_for_integration);

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
    ## Extract cell node coordinates.
    cell_node_coord_list = mesh_obj.mesh_nodes(mesh_obj.mesh_cells(e, :), :);
          
    ## Iterate over each quadrature point.
    for q = 1:quad_point_num
      ## The product of the Jacobian and the quadrature weights.
      JxW = qwts(q) * mesh_obj.cell_surface_metrics(e);
      ## Get the global coordinates for the current quadrature point
      ## and evaluate the functor at this point.
      global_coords = AreaToGlobalCoords(qpts(q, :), shape_functions_for_geometry, cell_node_coord_list);
      L1_norm += abs(functor(global_coords)) * JxW;
    endfor
  endfor

  L1_norm = L1_norm * triangle_unit_volume();
endfunction
