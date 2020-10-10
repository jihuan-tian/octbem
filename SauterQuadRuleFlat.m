function ret = SauterQuadRuleFlat(kernel_function,
				  kx_basis_function,
				  ky_basis_function,
				  kx_shape_functions_for_geometry,
				  ky_shape_functions_for_geometry,
				  kx_cell_index,
				  ky_cell_index,
				  panel_distance_matrix,
				  neighboring_type_matrix,
				  mesh_cells,
				  mesh_nodes,
				  nx, ny,
				  Jx_functor, Jy_functor,
				  mesh_size_estimate)
  ## SauterQuadRuleFlat - Perform integration on flat panels using the
  ## Sauter's method. The integral is
  ## \f[
  ## \int_{\hat{K} \times \hat{K}} H_{loc}(\xi, \eta) k(\xi, \eta) \intd s(\xi) \intd s(\eta),
  ## \f]
  ## 
  ## where \f$H_{loc}(\xi, \eta) = g_K(\xi) g_K(\eta)
  ## \hat{\varphi}(\xi) \hat{\psi}(\eta)\f$ and \f$k_{loc}(x, y)\f$ is
  ## the kernel function depending on global coordinates with \f$x, y
  ## \in \mathbb{R}^3\f$.
  ##
  ## Note: The integration of the product of \f$H(x, y)\f$ and \f$k(x,
  ## y)\f$ will be pulled back to the reference cell \f$\hat{K}\f$
  ## first, then to the master cell and will further be pulled back
  ## onto the parameter space \f$(\omega, \eta)\f$ for the numerical
  ## quadrature using Sauter's method.
  ##
  ## @param kernel_function The handle for the kernel function depending on
  ## global coordinates.
  ## @param kx_basis_function The test function associated with a supporting point
  ## in the reference cell. Note: the set of all supporting points is used for
  ## describing the function distribution instead of the real cell geometry.
  ## @param ky_basis_function The test function associated with a supporting point
  ## in the reference cell. Note: the set of all supporting points is used for
  ## describing the function distribution instead of the real cell geometry.
  ## @param kx_shape_functions_for_geometry A list of geometric shape functions
  ## for cell \f$K_x\f$.
  ## @param ky_shape_functions_for_geometry A list of geometric shape functions
  ## for cell \f$K_y\f$.
  ## @param kx_cell_index The cell index for \f$K_x\f$.
  ## @param ky_cell_index The cell index for \f$K_y\f$.
  ## @param panel_distance_matrix Distance between each pair of panels.
  ## @param neighboring_type_matrix Neighboring type between each pair of panels.
  ## @param mesh_cells All cells in the mesh.
  ## @param mesh_nodes All nodes in the mesh.
  ## @param nx The cell normal vector for \f$K_x\f$.
  ## @param ny The cell normal vector for \f$K_x\f$.
  ## @param Jx_functor The functor depending on area coordinates for
  ## calculating the Jacobian determinant for \f$K_x\f$.
  ## @param Jy_functor The functor depending on area coordinates for
  ## calculating the Jacobian determinant for \f$K_y\f$.
  ## @param mesh_size_estimate An estimate of the mesh size \f$h\f$, which is
  ## used for calculating quadrature order.
  ## @param ret The quadrature value.

  global sauter_same_panel_4d_qpts sauter_same_panel_4d_qwts sauter_common_edge_4d_qpts sauter_common_edge_4d_qwts sauter_common_vertex_4d_qpts sauter_common_vertex_4d_qwts sauter_regular_4d_qpts sauter_regular_4d_qwts;

  ## Determination of the quadrangle cell neighboring relationship.
  cell_neighboring_type = neighboring_type_matrix(kx_cell_index, ky_cell_index);

  ## Get the list of node indices stored in the two cells.
  kx_cell_node_indices = mesh_cells(kx_cell_index, :);
  ky_cell_node_indices = mesh_cells(ky_cell_index, :);

  ## Determine the polynomial order for describing the geometry of the two cells.
  kx_cell_order = QuadGeometryOrder(length(kx_cell_node_indices));
  ky_cell_order = QuadGeometryOrder(length(ky_cell_node_indices));

  ## Generate the matrix storing the standard ordered numberings of the node
  ## indices for the two cells in the standard order.
  kx_std_node_indices_matrix = GenQuadNodeIndicesMatrix(1:NumberOfNodesOnQuadrangle(kx_cell_order));
  ky_std_node_indices_matrix = GenQuadNodeIndicesMatrix(1:NumberOfNodesOnQuadrangle(ky_cell_order));
  
  ## Generate the matrix storing the node indices for the two cells in the
  ## standard order.
  kx_original_node_indices_matrix = GenQuadNodeIndicesMatrix(kx_cell_node_indices);
  ky_original_node_indices_matrix = GenQuadNodeIndicesMatrix(ky_cell_node_indices);
  
  ## Fixed quadrature order is adopted at the moment, so there is no
  ## quadrature order calculation like that in Erichsen quadrature
  ## rule.

  ret = 0;

  switch (cell_neighboring_type)
    case 1			# Identical
      ## fprintf(stderr(), "SauterQuadRule: same panel case!\n");
      
      ## Extract node coordinates.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices, :);
      
      ret = SauterQuadSamePanelFlat(kernel_function, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx_functor, Jy_functor, sauter_same_panel_4d_qpts, sauter_same_panel_4d_qwts);
    case 2			# Common edge
      ## fprintf(stderr(), "SauterQuadRule: common edge case!\n");
      
      ## Reverse the zigzag traversing order of the node indices for cell
      ## \f$K_y\f$.
      ky_cell_node_indices_reversed = GenPermutatedQuadNodeIndicesBackward(1, ky_original_node_indices_matrix);
      ## Recalculate the matrix for storing the reversed node indices.
      ky_reversed_node_indices_matrix = GenQuadNodeIndicesMatrix(ky_cell_node_indices_reversed);
      ## Get corner node indices for the two cells.
      kx_corner_node_indices = kx_cell_node_indices(GetQuadCornerNodeIndices(kx_cell_order));
      ky_corner_node_indices = ky_cell_node_indices_reversed(GetQuadCornerNodeIndices(ky_cell_order));
      ## Take the intersection of the corner node indices for the two cells.
      corner_node_index_intersection = intersect(kx_corner_node_indices, ky_corner_node_indices);
      if (length(corner_node_index_intersection) != 2)
	error("There should be two intersection nodes for the common edge case!");
      endif
      
      ## Get the ordering No. of the corner node in Kx from which the zigzag traversing should start.
      first_intersect_corner_node_idx_in_kx_corner_node_indices = find(kx_corner_node_indices == corner_node_index_intersection(1));
      second_intersect_corner_node_idx_in_kx_corner_node_indices = find(kx_corner_node_indices == corner_node_index_intersection(2));

      ## Check if the next corner node index belongs to the intersection in Kx.
      if ((mod(first_intersect_corner_node_idx_in_kx_corner_node_indices, 4) + 1) == second_intersect_corner_node_idx_in_kx_corner_node_indices)
	kx_corner_node_No_for_starting_zigzag = first_intersect_corner_node_idx_in_kx_corner_node_indices;
	starting_intersect_corner_node_index = corner_node_index_intersection(1);
      else
	kx_corner_node_No_for_starting_zigzag = second_intersect_corner_node_idx_in_kx_corner_node_indices;
	starting_intersect_corner_node_index = corner_node_index_intersection(2);
      endif

      ## Get the ordering No. of the corner node in Ky from which the zigzag traversing should start.
      ky_corner_node_No_for_starting_zigzag = find(ky_corner_node_indices == starting_intersect_corner_node_index);

      ## Calculate the permutation indices for the node indices in \f$K_x\f$.
      kx_node_permutation_indices = GenPermutatedQuadNodeIndicesForward(kx_corner_node_No_for_starting_zigzag, kx_std_node_indices_matrix);
      ## Permute the node indices for \f$K_x\f$ in the forward zigzag order.
      kx_cell_node_indices_perm = kx_cell_node_indices(kx_node_permutation_indices);
      ## Calculate the permutation indices for the reversed node indices in \f$K_y\f$.
      ky_node_permutation_indices = GenPermutatedQuadNodeIndicesForward(ky_corner_node_No_for_starting_zigzag, ky_std_node_indices_matrix);
      ## Permute the node indices for \f$K_y\f$ still in the forward zigzag
      ## order, because the node indices have already been reversed at the
      ## beginning.
      ky_cell_node_indices_perm = ky_cell_node_indices_reversed(ky_node_permutation_indices);
      
      ## Extract coordinates for ordered cell nodes.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices_perm, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices_perm, :);

      ## Overwrite the input Jx and Jy because the list of cell
      ## nodes are permutated.
      Jx_functor = @(kx_area_coord) GlobalSurfaceMetricOn3DQuad(kx_area_coord, kx_cell_node_coord_list);
      Jy_functor = @(ky_area_coord) GlobalSurfaceMetricOn3DQuad(ky_area_coord, ky_cell_node_coord_list);

      ## 2020-10-04: The following reordering of shape functions is redundant.
      ## ## It should be noted that the shape functions for describing the geometry
      ## ## of the cells should also be permuted.
      ## kx_shape_functions_for_geometry_perm = kx_shape_functions_for_geometry(kx_node_permutation_indices);
      ## ## Reverse the ordering of the shape functions on cell Ky.
      ## ky_shape_functions_for_geometry_reversed = ky_shape_functions_for_geometry(end:(-1):1);
      ## ky_shape_functions_for_geometry_perm = ky_shape_functions_for_geometry_reversed(ky_node_permutation_indices);
      ## ret = SauterQuadCommonEdgeFlat(kernel_function, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry_perm, ky_shape_functions_for_geometry_perm, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx_functor, Jy_functor, sauter_common_edge_4d_qpts, sauter_common_edge_4d_qwts);

      ret = SauterQuadCommonEdgeFlat(kernel_function, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx_functor, Jy_functor, sauter_common_edge_4d_qpts, sauter_common_edge_4d_qwts);
    case 3			# Common vertex
      ## fprintf(stderr(), "SauterQuadRule: common vertex case!\n");
      
      ## Get corner node indices for the two cells.
      kx_corner_node_indices = kx_cell_node_indices(GetQuadCornerNodeIndices(kx_cell_order));
      ky_corner_node_indices = ky_cell_node_indices(GetQuadCornerNodeIndices(ky_cell_order));
      
      ## Take the intersection of the two cells of corner node indices.
      corner_node_index_intersection = intersect(kx_corner_node_indices, ky_corner_node_indices);
      if (length(corner_node_index_intersection) != 1)
	error("There should be only one intersection node for the common vertex case!");
      endif

      ## The common vertex is the starting corner point for zigzag traversing.
      kx_corner_node_No_for_starting_zigzag = find(kx_corner_node_indices == corner_node_index_intersection);
      ky_corner_node_No_for_starting_zigzag = find(ky_corner_node_indices == corner_node_index_intersection);

      ## Calculate the permutation indices for the node indices in \f$K_x\f$.
      kx_node_permutation_indices = GenPermutatedQuadNodeIndicesForward(kx_corner_node_No_for_starting_zigzag, kx_std_node_indices_matrix);
      ## Permutate the node indices for \f$K_x\f$ in the forward zigzag order.
      kx_cell_node_indices_perm = kx_cell_node_indices(kx_node_permutation_indices);
      ## Calculate the permutation indices for the node indices in \f$K_y\f$.
      ky_node_permutation_indices = GenPermutatedQuadNodeIndicesForward(ky_corner_node_No_for_starting_zigzag, ky_std_node_indices_matrix);
      ## Permutate the node indices for \f$K_y\f$ in the forward zigzag order.
      ky_cell_node_indices_perm = ky_cell_node_indices(ky_node_permutation_indices);

      ## Extract coordinates for ordered cell nodes.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices_perm, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices_perm, :);

      ## Overwrite the input Jx and Jy because the list of cell
      ## nodes are permutated.
      Jx_functor = @(kx_area_coord) GlobalSurfaceMetricOn3DQuad(kx_area_coord, kx_cell_node_coord_list);
      Jy_functor = @(ky_area_coord) GlobalSurfaceMetricOn3DQuad(ky_area_coord, ky_cell_node_coord_list);

      ## 2020-10-04: The following reordering of shape functions is redundant.
      ## It should be noted that the shape functions for describing the geometry
      ## of the cells should also be permuted.
      ## kx_shape_functions_for_geometry_perm = kx_shape_functions_for_geometry(kx_node_permutation_indices);
      ## ky_shape_functions_for_geometry_perm = ky_shape_functions_for_geometry(ky_node_permutation_indices);
      ## ret = SauterQuadCommonVertexFlat(kernel_function, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry_perm, ky_shape_functions_for_geometry_perm, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx_functor, Jy_functor, sauter_common_vertex_4d_qpts, sauter_common_vertex_4d_qwts);

      ret = SauterQuadCommonVertexFlat(kernel_function, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx_functor, Jy_functor, sauter_common_vertex_4d_qpts, sauter_common_vertex_4d_qwts);
    case 4			# Regular
      ## fprintf(stderr(), "SauterQuadRule: regular case!\n");
      
      ## Extract coordinates for ordered cell nodes.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices, :);
      
      ## Get the distance between the two panels.
      panel_distance = panel_distance_matrix(kx_cell_index, ky_cell_index);
            
      ret = SauterQuadRegularFlat(kernel_function, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx_functor, Jy_functor, sauter_regular_4d_qpts, sauter_regular_4d_qwts);
  endswitch
endfunction
