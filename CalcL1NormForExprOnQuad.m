function L1_norm = CalcL1NormForExprOnQuad(functor, shape_functions_for_geometry,
					   cell_indices_for_integration, mesh_obj, qpts, qwts)
  ## CalcL1NormForExprOnQuad - Calculate the L1-norm of the given analytical
  ## expression for quadrangle mesh.
  ## @param functor The functor defining the analytical expression
  ## dependent on the global coordinates.
  ## @param shape_functions_for_geometry The shape function space for
  ## describing the cell geometry.
  ## @param cell_indices_for_integration The list of cell indices for
  ## carrying out the integration.
  ## @param mesh_obj The object holding all the mesh data.
  ## @param qpts Quadrature points on the quadrangle.
  ## @param qwts Quadrature weights on the quadrature.
  ## @param L1_norm Calculated L1 norm.

  number_of_cells_for_integration = length(cell_indices_for_integration);

  quad_point_num = size(qwts, 1);

  L1_norm = 0;

  ## Iterate over each cell for the integration.
  for m = 1:number_of_cells_for_integration
    ## Current cell index
    e = cell_indices_for_integration(m);
    ## Extract cell node coordinates.
    cell_node_coord_list = mesh_obj.mesh_nodes(mesh_obj.mesh_cells(e, :), :);

    ## Construct surface metric functor on the current cell.
    J = @(area_coord) GlobalSurfaceMetricOn3DQuad(area_coord, cell_node_coord_list);
          
    ## Iterate over each quadrature point.
    for q = 1:quad_point_num
      ## The product of the Jacobian and the quadrature weights.
      JxW = qwts(q) * J(qpts(q, :));
      ## Get the global coordinates for the current quadrature point
      ## and evaluate the functor at this point.
      global_coords = AreaToGlobalCoords(qpts(q, :), shape_functions_for_geometry, cell_node_coord_list);
      L1_norm += abs(functor(global_coords)) * JxW;
    endfor
  endfor
endfunction
