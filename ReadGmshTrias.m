function gmsh_mesh_obj = ReadGmshTrias(gmsh_msh_file_name)
  ## ReadGmshTrias - Read the triangle cells from Gmsh mesh file.
  ## @param gmsh_msh_file_name The file name of the Gmsh msh file.
  ## @param gmsh_mesh_obj The struct object containg Gmsh mesh data, which
  ## includes the following fields:
  ##   mesh_nodes: A list of nodes stored in N*M matrix, with N as the number of
  ##   nodes, M as the space dimension.
  ##   mesh_cells: A list of triangle cells stored in N*M matrix, with N as the
  ##   number of cells, M as the number of nodes in a cell.
  ##   number_of_nodes: Total number of nodes.
  ##   space_dim: The space dimension used to define the node coordinates.
  ##   number_of_cells: Total number of cells in the mesh.
  ##   number_of_nodes_in_cells: Total number of nodes in each cell, which
  ##   determines the polynomial order of the shape functions for describing the
  ##   cell geometry.
  ##   cell_geom_orders: The polynomial orders of the shape functions for
  ##   describing the geometry of all cells.
  ##   corner_node_local_indices_in_cells: Local indices of the corner nodes in
  ##   all cells, where it is assumed that all the nodes are in the forward
  ##   zigzag order.
  ##   cell_ranges: The ranges of all cells.
  ##   min_cell_range: The minimum mesh size \f$h\f$.
  ##   cell_surface_metrics: Calculate the surface metric from the reference
  ##   cell to each real cell, which is the squared sum of the three
  ##   two-coordinate-component Jacobians.
  ##   cell_normal_vectors: Calculate the surface normal vector for each cell.

  ## Get the mesh data size.
  [gmsh_mesh_obj.number_of_nodes, gmsh_mesh_obj.space_dim, gmsh_mesh_obj.number_of_cells, number_of_nodes_in_cell] = gmsh_size_read(gmsh_msh_file_name);
  gmsh_mesh_obj.number_of_nodes_in_cells = ones(gmsh_mesh_obj.number_of_cells, 1) * number_of_nodes_in_cell;
  ## Read mesh data from the file.
  [gmsh_mesh_obj.mesh_nodes, gmsh_mesh_obj.mesh_cells] = gmsh_data_read(gmsh_msh_file_name, gmsh_mesh_obj.space_dim, gmsh_mesh_obj.number_of_nodes, number_of_nodes_in_cell, gmsh_mesh_obj.number_of_cells);
  ## Calculate the shape function orders for describing the geometry of all cells.
  gmsh_mesh_obj.cell_geom_orders = ones(gmsh_mesh_obj.number_of_cells, 1) * TriaGeometryOrder(number_of_nodes_in_cell);
  ## Local indices of the corner nodes in all cells, where it is assumed that all
  ## the nodes are in the forward zigzag order.
  gmsh_mesh_obj.corner_node_local_indices_in_cells = repmat(GetTriaCornerNodeIndices(gmsh_mesh_obj.cell_geom_orders(1)), gmsh_mesh_obj.number_of_cells, 1);
  
  gmsh_mesh_obj.cell_ranges = zeros(gmsh_mesh_obj.number_of_cells, 1);
  gmsh_mesh_obj.cell_surface_metrics = zeros(gmsh_mesh_obj.number_of_cells, 1);
  gmsh_mesh_obj.cell_normal_vectors = zeros(gmsh_mesh_obj.number_of_cells, gmsh_mesh_obj.space_dim);
  gmsh_mesh_obj.cell_circumcenters = zeros(gmsh_mesh_obj.number_of_cells, gmsh_mesh_obj.space_dim);
  gmsh_mesh_obj.cell_cogs = zeros(gmsh_mesh_obj.number_of_cells, gmsh_mesh_obj.space_dim);

  for m = 1:gmsh_mesh_obj.number_of_cells
    ## Get the corner node coordinates for the current cell.
    cell_corner_nodes = gmsh_mesh_obj.mesh_nodes(gmsh_mesh_obj.mesh_cells(m, gmsh_mesh_obj.corner_node_local_indices_in_cells(m, :)), :);
    ## Estimate the triangle range \f$h\f$ for each cell.
    gmsh_mesh_obj.cell_ranges(m) = Tria3DRange(cell_corner_nodes);
    ## Calculate the surface metric from the reference cell to each real cell.
    gmsh_mesh_obj.cell_surface_metrics(m) = GlobalSurfaceMetricOn3DFlatTria(cell_corner_nodes);
    ## Calculate the normal vector of each cell.
    gmsh_mesh_obj.cell_normal_vectors(m, :) = SurfaceNormalOn3DFlatTria(cell_corner_nodes);
    ## Calculate the circumcenter of the current cell.
    switch (gmsh_mesh_obj.space_dim)
      case 2
	gmsh_mesh_obj.cell_circumcenters(m, :) = Tria2DCircumcenter(cell_corner_nodes);
      case 3
	gmsh_mesh_obj.cell_circumcenters(m, :) = Tria3DCircumcenter(cell_corner_nodes);
    endswitch

    ## Calculate the center of gravity of the current cell.
    gmsh_mesh_obj.cell_cogs(m, :) = TriaCenterOfGravity(cell_corner_nodes);
  endfor
  
  gmsh_mesh_obj.min_cell_range = min(gmsh_mesh_obj.cell_ranges);
endfunction
