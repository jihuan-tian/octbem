function dist_matrix = CalcPanelDistanceMatrix(mesh_data)
  ## CalcPanelDistanceMatrix - Calculate the distance between each
  ## pair of panels.
  ## @param mesh_data The object holding mesh data.
  ## @param dist_matrix Distance between each pair of panels.

  dist_matrix = zeros(mesh_data.number_of_cells, mesh_data.number_of_cells);

  for m = 1:mesh_data.number_of_cells
    kx_cell_node_indices = mesh_data.mesh_cells(m, :);
    kx_cell_node_coord_list = mesh_data.mesh_nodes(kx_cell_node_indices, :);
    
    kx_cell_order = TriaGeometryOrder(length(kx_cell_node_indices));
    kx_corner_node_indices = kx_cell_node_indices(GetTriaCornerNodeIndices(kx_cell_order));
    kx_corner_node_coord_list = mesh_data.mesh_nodes(kx_corner_node_indices, :);

    for n = (m+1):mesh_data.number_of_cells
      ky_cell_node_indices = mesh_data.mesh_cells(n, :);
      ky_cell_node_coord_list = mesh_data.mesh_nodes(ky_cell_node_indices, :);

      ky_cell_order = TriaGeometryOrder(length(ky_cell_node_indices));
      ky_corner_node_indices = ky_cell_node_indices(GetTriaCornerNodeIndices(ky_cell_order));
      ky_corner_node_coord_list = mesh_data.mesh_nodes(ky_corner_node_indices, :);

      dist_matrix(m, n) = DistOfConvexPolygons(kx_cell_node_coord_list, kx_corner_node_coord_list, ky_cell_node_coord_list, ky_corner_node_coord_list);
      dist_matrix(n, m) = dist_matrix(m, n);
    endfor
  endfor
endfunction
