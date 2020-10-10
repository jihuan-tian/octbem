function ret = ErichsenQuadRuleFlat(kernel_function, kernel_singularity_order,
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
				    Jx, Jy,
				    mesh_size_estimate,
				    sobolev_function_space_order,
				    basis_function_polynomial_order,
				    galerkin_estimate_norm_index)
  ## ErichsenQuadRuleFlat - Perform integration on flat panels using the method
  ## in Erichsen1996Efficient. The integral is
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
  ## quadrature using Erichsen's method.
  ##
  ## @param kernel_function The handle for the kernel function depending on
  ## global coordinates.
  ## @param kernel_singularity_order Singularity order of the kernel to be
  ## integrated.
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
  ## @param Jx The surface metric (Jacobian) for \f$K_x\f$ from reference to real cell.
  ## @param Jy The surface metric (Jacobian) for \f$K_y\f$ from reference to real cell.
  ## @param mesh_size_estimate An estimate of the mesh size \f$h\f$, which is
  ## used for calculating quadrature order.
  ## @param sobolev_function_space_order The Sobolev space order for the
  ## function space, i.e. \f$s_1\f$
  ## @param basis_function_polynomial_order The polynomial order of the adopted
  ## basis function.
  ## @param galerkin_estimate_norm_index The index of the Sobolev norm, which is
  ## used for measuring the error in Galerkin esitmate. It is equal to \f$s_1 -
  ## t\f$. When it is zero, the error is \f$L^2\f$.
  ## @param ret The quadrature value.

  ## Determination of the triangular cell neighboring relationship.
  cell_neighboring_type = neighboring_type_matrix(kx_cell_index, ky_cell_index);

  ## Get the list of node indices stored in the two cells.
  kx_cell_node_indices = mesh_cells(kx_cell_index, :);
  ky_cell_node_indices = mesh_cells(ky_cell_index, :);

  ## Determine the polynomial order for describing the geometry of the two cells.
  kx_cell_order = TriaGeometryOrder(length(kx_cell_node_indices));
  ky_cell_order = TriaGeometryOrder(length(ky_cell_node_indices));

  ## Generate the matrix storing the standard ordered numberings of the node
  ## indices for the two cells in the standard order.
  kx_std_node_indices_matrix = GenTriaNodeIndicesMatrix(1:NumberOfNodesOnTriangle(kx_cell_order));
  ky_std_node_indices_matrix = GenTriaNodeIndicesMatrix(1:NumberOfNodesOnTriangle(ky_cell_order));
  
  ## Generate the matrix storing the node indices for the two cells in the
  ## standard order.
  kx_original_node_indices_matrix = GenTriaNodeIndicesMatrix(kx_cell_node_indices);
  ky_original_node_indices_matrix = GenTriaNodeIndicesMatrix(ky_cell_node_indices);

  ## Calculate the quadrature order, actually, the number of quadrature
  ## points.
  [norder_for_eta, norder_for_omega] = ErichsenSingularQuadOrder(basis_function_polynomial_order, sobolev_function_space_order, mesh_size_estimate);

  ret = 0;

  switch (cell_neighboring_type)
    case 1			# Identical
      ## fprintf(stderr(), "Erichsen1996Efficient: same panel case!\n");
      
      ## Extract node coordinates.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices, :);
      
      ret = ErichsenQuadSamePanelFlat(kernel_function, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx, Jy);
    case 2			# Common edge
      ## fprintf(stderr(), "Erichsen1996Efficient: common edge case!\n");
      
      ## Reverse the zigzag traversing order of the node indices for cell
      ## \f$K_y\f$.
      ky_cell_node_indices_reversed = GenPermutatedTriaNodeIndicesBackward(1, ky_original_node_indices_matrix);
      ## Recalculate the lower triangular matrix for storing the reversed node indices.
      ky_reversed_node_indices_matrix = GenTriaNodeIndicesMatrix(ky_cell_node_indices_reversed);
      ## Get corner node indices for the two cells.
      kx_corner_node_indices = kx_cell_node_indices(GetTriaCornerNodeIndices(kx_cell_order));
      ky_corner_node_indices = ky_cell_node_indices_reversed(GetTriaCornerNodeIndices(ky_cell_order));
      ## Take the intersection of the corner node indices for the two cells.
      corner_node_index_intersection = intersect(kx_corner_node_indices, ky_corner_node_indices);
      if (length(corner_node_index_intersection) != 2)
	error("There should be two intersection nodes for the common edge case!");
      endif
      
      ## Get the ordering No. of the corner node in Kx which does not lie on the common edge.
      kx_not_shared_corner_node_No = find((kx_corner_node_indices != corner_node_index_intersection(1)) & (kx_corner_node_indices != corner_node_index_intersection(2)));
      ## Get the ordering No. of the corner node in Kx from which the zigzag traversing should start.
      if (kx_not_shared_corner_node_No != 1)
	kx_corner_node_No_for_starting_zigzag = kx_not_shared_corner_node_No - 1;
      else
	kx_corner_node_No_for_starting_zigzag = 3;
      endif

      ## Get the ordering No. of the corner node in Ky which does not lie on the common edge.
      ky_not_shared_corner_node_No = find((ky_corner_node_indices != corner_node_index_intersection(1)) & (ky_corner_node_indices != corner_node_index_intersection(2)));
      ## Get the ordering No. of the corner node in Ky from which the zigzag traversing should start.
      if (ky_not_shared_corner_node_No != 1)
	ky_corner_node_No_for_starting_zigzag = ky_not_shared_corner_node_No - 1;
      else
	ky_corner_node_No_for_starting_zigzag = 3;
      endif

      ## Calculate the permutation indices for the node indices in \f$K_x\f$.
      kx_node_permutation_indices = GenPermutatedTriaNodeIndicesForward(kx_corner_node_No_for_starting_zigzag, kx_std_node_indices_matrix);
      ## Permute the node indices for \f$K_x\f$ in the forward zigzag order.
      kx_cell_node_indices_perm = kx_cell_node_indices(kx_node_permutation_indices);
      ## Calculate the permutation indices for the reversed node indices in \f$K_y\f$.
      ky_node_permutation_indices = GenPermutatedTriaNodeIndicesForward(ky_corner_node_No_for_starting_zigzag, ky_std_node_indices_matrix);
      ## Permute the node indices for \f$K_y\f$ still in the forward zigzag
      ## order, because the node indices have already been reversed at the
      ## beginning.
      ky_cell_node_indices_perm = ky_cell_node_indices_reversed(ky_node_permutation_indices);
      
      ## Extract coordinates for ordered cell nodes.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices_perm, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices_perm, :);

      ## 2020-10-04: The following reordering of shape functions is redundant.
      ## ## It should be noted that the shape functions for describing the geometry
      ## ## of the cells should also be permuted.
      ## kx_shape_functions_for_geometry_perm = kx_shape_functions_for_geometry(kx_node_permutation_indices);
      ## ## Reverse the ordering of the shape functions on cell Ky.
      ## ky_shape_functions_for_geometry_reversed = ky_shape_functions_for_geometry(end:(-1):1);
      ## ky_shape_functions_for_geometry_perm = ky_shape_functions_for_geometry_reversed(ky_node_permutation_indices);
      ## ret = ErichsenQuadCommonEdgeFlat(kernel_function, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry_perm, ky_shape_functions_for_geometry_perm, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx, Jy);

      ret = ErichsenQuadCommonEdgeFlat(kernel_function, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx, Jy);
    case 3			# Common vertex
      ## fprintf(stderr(), "Erichsen1996Efficient: common vertex case!\n");
      
      ## Get corner node indices for the two cells.
      kx_corner_node_indices = kx_cell_node_indices(GetTriaCornerNodeIndices(kx_cell_order));
      ky_corner_node_indices = ky_cell_node_indices(GetTriaCornerNodeIndices(ky_cell_order));
      
      ## Take the intersection of the two cells of corner node indices.
      corner_node_index_intersection = intersect(kx_corner_node_indices, ky_corner_node_indices);
      if (length(corner_node_index_intersection) != 1)
	error("There should be only one intersection node for the common vertex case!");
      endif

      kx_common_corner_node_No = find(kx_corner_node_indices == corner_node_index_intersection);
      ky_common_corner_node_No = find(ky_corner_node_indices == corner_node_index_intersection);
      kx_corner_node_No_for_starting_zigzag = mod((kx_common_corner_node_No - 1) + 1, 3) + 1;
      ky_corner_node_No_for_starting_zigzag = mod((ky_common_corner_node_No - 1) + 1, 3) + 1;

      ## Calculate the permutation indices for the node indices in \f$K_x\f$.
      kx_node_permutation_indices = GenPermutatedTriaNodeIndicesForward(kx_corner_node_No_for_starting_zigzag, kx_std_node_indices_matrix);
      ## Permutate the node indices for \f$K_x\f$ in the forward zigzag order.
      kx_cell_node_indices_perm = kx_cell_node_indices(kx_node_permutation_indices);
      ## Calculate the permutation indices for the node indices in \f$K_y\f$.
      ky_node_permutation_indices = GenPermutatedTriaNodeIndicesForward(ky_corner_node_No_for_starting_zigzag, ky_std_node_indices_matrix);
      ## Permutate the node indices for \f$K_y\f$ in the forward zigzag order.
      ky_cell_node_indices_perm = ky_cell_node_indices(ky_node_permutation_indices);

      ## Extract coordinates for ordered cell nodes.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices_perm, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices_perm, :);

      ## 2020-10-04: The following reordering of shape functions is redundant.
      ## ## It should be noted that the shape functions for describing the geometry
      ## ## of the cells should also be permuted.
      ## kx_shape_functions_for_geometry_perm = kx_shape_functions_for_geometry(kx_node_permutation_indices);
      ## ky_shape_functions_for_geometry_perm = ky_shape_functions_for_geometry(ky_node_permutation_indices);
      ## ret = ErichsenQuadCommonVertexFlat(kernel_function, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry_perm, ky_shape_functions_for_geometry_perm, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx, Jy);
      
      ret = ErichsenQuadCommonVertexFlat(kernel_function, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx, Jy);
    case 4			# Regular
      ## fprintf(stderr(), "Erichsen1996Efficient: regular case!\n");
      
      ## Extract coordinates for ordered cell nodes.
      kx_cell_node_coord_list = mesh_nodes(kx_cell_node_indices, :);
      ky_cell_node_coord_list = mesh_nodes(ky_cell_node_indices, :);
      
      ## Get the distance between the two panels.
      panel_distance = panel_distance_matrix(kx_cell_index, ky_cell_index);
      
      ## Calculate the number of quadrature points.
      norder = ErichsenRegularQuadOrder(kernel_singularity_order, basis_function_polynomial_order, sobolev_function_space_order, mesh_size_estimate, panel_distance, galerkin_estimate_norm_index);
      
      ret = ErichsenQuadRegularFlat(kernel_function, norder, kx_basis_function, ky_basis_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, nx, ny, Jx, Jy);
  endswitch
endfunction
