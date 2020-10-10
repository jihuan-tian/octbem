function ret = SauterQuadCommonEdgeFlat(kernel_function, 
					kx_basis_function,
					ky_basis_function,
					kx_shape_functions_for_geometry,
					ky_shape_functions_for_geometry,
					kx_cell_node_coord_list,
					ky_cell_node_coord_list,
					nx, ny,
					Jx_functor, Jy_functor,
					qpts, qwts)
  ## SauterQuadCommonEdgeFlat - Perform the common edge integration
  ## using the method in Sauter's BEM, where the quadrangle is flat.
  ## The integral is
  ## \f[
  ## \int_{\hat{K} \times \hat{K}} H_{loc}(\xi, \eta) k_{loc}(\xi, \eta) \intd s(\xi) \intd s(\eta),
  ## \f]
  ## 
  ## where \f$H_{loc}(\xi, \eta) = g_K(\xi) g_K(eta)
  ## \hat{\varphi}(\xi) \hat{\psi}(\eta)\f$ and \f$k_{loc}(x, y)\f$ is
  ## the kernel function depending on global coordinates with \f$x, y
  ## \in \mathbb{R}^3\f$.
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
  ## @param kx_cell_node_coord_list A list of node global coordinates for the
  ## cell \f$K_x\f$.
  ## @param ky_cell_node_coord_list A list of node global coordinates for the
  ## cell \f$K_y\f$.
  ## @param nx The cell normal vector for \f$K_x\f$.
  ## @param ny The cell normal vector for \f$K_y\f$.
  ## @param Jx_functor The functor depending on area coordinates for
  ## calculating the Jacobian determinant for \f$K_x\f$.
  ## @param Jy_functor The functor depending on area coordinates for
  ## calculating the Jacobian determinant for \f$K_y\f$.
  ## @param qpts Quadrature points for the 4d parameter space f$(0,1)^4\f$.
  ## @param qwts Quadrature weights for the 4d parameter space f$(0,1)^4\f$.
  ## @param ret The quadrature value.
  
  ## Generate the pullback of the kernel function \f$k(x, y)\f$ to the reference
  ## cell \f$\hat{K}_x \times \hat{K}_y\f$.
  k_loc = @(kx_area_coord, ky_area_coord) BEMKernelPullbackFlat(kx_area_coord, ky_area_coord, nx, ny, kernel_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list);
  
  ## Generate the pullback of \f$H_loc\f$ to the reference cell
  ## \f$\hat{K}_x \times \hat{K}_y\f$.
  H_loc = @(kx_area_coord, ky_area_coord) Jx_functor(kx_area_coord) * Jy_functor(ky_area_coord) * kx_basis_function(kx_area_coord) * ky_basis_function(ky_area_coord);
  
  ## Generate the integrand which combines both \f$H_{loc}\f$ and \f$K_{loc}\f$
  ## on the reference cells.
  k3 = @(kx_area_coord, ky_area_coord) H_loc(kx_area_coord, ky_area_coord) * k_loc(kx_area_coord, ky_area_coord);
  
  ## Pullback the integrand \f$k_3\f$ to the 4d parameter domain. The
  ## parameter space is \f$[\xi, \eta_1, \eta_2, \eta_3]\f$.
  k3_tilde = @(p) p(:,1).^2 .* (1 - p(:,1)) .* (
	       k3([(1 - p(:,1)) .* p(:,4) + p(:,1), p(:,1) .* p(:,3)], ...
		  [(1 - p(:,1)) .* p(:,4), p(:,1) .* p(:,2)]) + ...
	       k3([(1 - p(:,1)) .* p(:,4), p(:,1) .* p(:,3)], ...
		  [p(:,1) + (1 - p(:,1)) .* p(:,4), p(:,1) .* p(:,2)])) + ...
	     p(:,1).^2 .* (1 - p(:,1) .* p(:,2)) .* (
	       k3([(1 - p(:,1) .* p(:,2)) .* p(:,4) + p(:,1) .* p(:,2), ...
		   p(:,1) .* p(:,3)], 
		  [(1 - p(:,1) .* p(:,2)) .* p(:,4), p(:,1)]) + 
	       k3([(1 - p(:,1) .* p(:,2)) .* p(:,4) + p(:,1) .* p(:,2), p(:,1)],
		  [(1 - p(:,1) .* p(:,2)) .* p(:,4), p(:,1) .* p(:,3)]) +
	       k3([(1 - p(:,1) .* p(:,2)) .* p(:,4), p(:,1) .* p(:,3)],
		  [(1 - p(:,1) .* p(:,2)) .* p(:,4) + p(:,1) .* p(:,2), p(:,1)]) +
	       k3([(1 - p(:,1) .* p(:,2)) .* p(:,4), p(:,1)],
		  [(1 - p(:,1) .* p(:,2)) .* p(:,4) + p(:,1) .* p(:,2), ...
		   p(:,1) .* p(:,3)]));
  
  ## Apply 4d Gauss-Legendre quadrature rule to the integrand.
  ret = ApplyQuadRule(k3_tilde, qpts, qwts);
endfunction
