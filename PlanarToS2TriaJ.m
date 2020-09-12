function jacobian_mats = PlanarToS2TriaJ(area_coordinates, shape_functions_for_geometry, cell_node_coord_list, sphere_center, sphere_radius)
  ## PlanarToS2TriaJ - Calculate the Jacobian matrix at the area coordinate.

  ## @param area_coordinates A list of area coordinates stored in a matrix with
  ## dimension N*M, where N is the number of area coordiantes, M is the number
  ## of coordinate components.
  ## @param shape_functions_for_geometry An array of shape function handles used for
  ## describing cell geometry, which are associated with supporting nodes in the
  ## cell.
  ## @param cell_node_coord_list A list of global coordinates for cell nodes with
  ## dimension N*M, where N is the number of cell nodes, M is the space
  ## dimension.
  ## @param sphere_center Center coordinates of the sphere manifold,
  ## which is stored in a row vector.
  ## @param sphere_radius The radius of the sphere model.
  ## @param jacobian_mats A cell array of Jacobian matrices at all points.

  number_of_points = size(area_coordinates, 1);
  jacobian_mats = cell(number_of_points, 1);
  global_coords = AreaToGlobalCoords(area_coordinates, shape_functions_for_geometry, cell_node_coord_list);

  for m = 1:number_of_points
    jacobian_mats{m} = zeros(3, 3);

    r = norm(global_coords(m, :) - sphere_center);
    r_cubic = r^3;

    ## Generate the Jacobian matrix.
    for i = 1:3
      for j = 1:i
	jacobian_mats{m}(i, j) = -(global_coords(m, i) - sphere_center(i)) * (global_coords(m, j) - sphere_center(j)) / r_cubic;
	jacobian_mats{m}(j, i) = jacobian_mats{m}(i, j);
      endfor
    endfor

    ## Modify the diagonal of the Jacobian matrix.
    for i = 1:3
      jacobian_mats{m}(i, i) += 1.0 / r;
    endfor

    ## Scale the Jacobian matrix
    jacobian_mats{m} *= sphere_radius;
  endfor
endfunction
