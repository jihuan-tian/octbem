function h = quadplot3d(mesh_cells, mesh_nodes, varargin)
  ## quadplot3d - Plot a list of quadrangles embedded in 3D space.
  ## @param mesh_cells A list of quadrangular cells, each of which is stored as an
  ## \f$N*M\f$ matrix. Each row of the matrix represents a cell, which is
  ## comprised of an ordered list of node indices.
  ## @param x An array of x components for the node coordinates.
  ## @param y An array of y components for the node coordinates.
  ## @param z An array of z components for the node coordinates.
  ## @param varargin Specification of line styles.
  ## @param h Figure handle.

  if (nargin < 2)
    print_usage();
  endif

  [number_of_cells, number_of_nodes_in_cell] = size(mesh_cells);
  [number_of_nodes, space_dim] = size(mesh_nodes);

  if (space_dim != 3)
    error("The space dimension should be 3!");
  endif
  
  ## Get the shape function order for describing the cell geometry.
  cell_geom_order = QuadGeometryOrder(number_of_nodes_in_cell);

  ## Get the corner node indices of the triangular cells.
  mesh_cells_corner_only = mesh_cells(:, GetQuadCornerNodeIndices(cell_geom_order));
  ## Generate the local node indices to form a loop of the quadrangle corner
  ## nodes.
  corner_node_indices_loop = mesh_cells_corner_only(:, [1, 2, 3, 4, 1]);
  
  plot3(mesh_nodes(:, 1)(corner_node_indices_loop)', mesh_nodes(:, 2)(corner_node_indices_loop)', mesh_nodes(:, 3)(corner_node_indices_loop)', varargin{:});

  if (nargout > 0)
    h = handle;
  endif
endfunction
