function ret = ErichsenQuadRegularFlat(kernel_function, norder, kx_basis_function,
				       ky_basis_function,
				       kx_shape_functions_for_geometry,
				       ky_shape_functions_for_geometry,
				       kx_cell_node_coord_list,
				       ky_cell_node_coord_list,
				       nx, ny,
				       Jx, Jy)
  ## ErichsenQuadRegularFlat - Perform the regular integration using the method in
  ## Erichsen1996Efficient, where the triangle is flat. The integral is
  ## \f[
  ## \int_{\hat{K} \times \hat{K}} H_{loc}(\xi, \eta) k(\xi, \eta) \intd s(\xi) \intd s(\eta),
  ## \f]
  ## 
  ## where \f$H_{loc}(\xi, \eta) = g_K(\xi) g_K(eta)
  ## \hat{\varphi}(\xi) \hat{\psi}(\eta)\f$ and \f$k_{loc}(x, y)\f$ is
  ## the kernel function depending on global coordinates with \f$x, y
  ## \in \mathbb{R}^3\f$.
  ##
  ## Note: The integration of the product of \f$H(x, y)\f$ and \f$k(x, y)\f$
  ## will be pulled back to the reference cell \f$\hat{K}\f$.
  ##
  ## @param kernel_function The handle for the kernel function depending on
  ## global coordinates.
  ## @param norder Number of quadrature points.
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
  ## @param kx_cell_node_coord_list A list of node global coordinates for the
  ## cell \f$K_x\f$.
  ## @param ky_cell_node_coord_list A list of node global coordinates for the
  ## cell \f$K_y\f$.
  ## @param nx The cell normal vector for \f$K_x\f$.
  ## @param ny The cell normal vector for \f$K_x\f$.
  ## @param Jx The surface metric (Jacobian) for \f$K_x\f$ from reference to real cell.
  ## @param Jy The surface metric (Jacobian) for \f$K_y\f$ from reference to real cell.
  ## @param ret The quadrature value.

  ## Generate the pullback of the kernel function \f$k(x, y)\f$ to the reference
  ## cell \f$\hat{K}_x \times \hat{K}_y\f$.
  k_loc = @(kx_area_coord, ky_area_coord) BEMKernelPullbackFlat(kx_area_coord, ky_area_coord, nx, ny, kernel_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list);
  
  ## Generate the pullback of \f$H_loc\f$ to the reference cell
  ## \f$\hat{K}_x \times \hat{K}_y\f$.
  H_loc = @(kx_area_coord, ky_area_coord) Jx * Jy * kx_basis_function(kx_area_coord) * ky_basis_function(ky_area_coord);
  
  ## Generate the integrand which combines both \f$H_{loc}\f$ and \f$K_{loc}\f$
  ## on the reference cells.
  full_integrand_on_Khat = @(kx_area_coord, ky_area_coord) H_loc(kx_area_coord, ky_area_coord) * k_loc(kx_area_coord, ky_area_coord);

  ## #############################################################################
  ## Integrate the integrand using the normal Gauss quadrature on the product
  ## space of two simplices.
  ## #############################################################################
  ## Generate the quadrature rule on unit triangle (2D simplex).
  [qpts1, qpts2, qwts] = triangle_unit_product_set(norder);

  ## Integrate on \f$K_y\f$.
  integral_on_ky = @(kx_area_coord) ApplyQuadRule(@(ky_area_coord) full_integrand_on_Khat(kx_area_coord, ky_area_coord), [qpts1', qpts2'], qwts') * triangle_unit_volume();

  ## Integrate on \f$K_x\f$.
  ret = ApplyQuadRule(integral_on_ky, [qpts1', qpts2'], qwts') * triangle_unit_volume();
endfunction
