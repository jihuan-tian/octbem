function ret = ErichsenQuadSamePanelS2(kernel_function, norder_for_eta,
				       norder_for_omega, kx_basis_function,
				       ky_basis_function,
				       kx_shape_functions_for_geometry,
				       ky_shape_functions_for_geometry,
				       kx_cell_node_coord_list,
				       ky_cell_node_coord_list,
				       nx_functor, ny_functor,
				       Jx_functor, Jy_functor,
				       sphere_center, sphere_radius)
  ## ErichsenQuadSamePanelS2 - Perform the same panel integration using the method in
  ## Erichsen1996Efficient. The integral is
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
  ## will be pulled back to the reference cell \f$\hat{K}\f$ first and will
  ## further be pulled back onto the parameter space \f$(\omega, \eta)\f$ for
  ## the numerical quadrature using Erichsen's method.
  ##
  ## @param kernel_function The handle for the kernel function depending on
  ## global coordinates.
  ## @param norder_for_eta The number of quadrature points for \f$eta\f$ on \f$R_{d_1}\f$.
  ## @param norder_for_omega The number of quadrature points for \f$\omega\f$ on \f$T_{d_2}\f$.
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
  ## @param nx_functor The functor depending on area coordinates for
  ## calculating the cell normal vectors for \f$K_x\f$.
  ## @param ny_functor The functor depending on area coordinates for
  ## calculating the cell normal vectors for \f$K_y\f$.
  ## @param Jx_functor The functor depending on area coordinates for
  ## calculating the Jacobian determinant for \f$K_x\f$.
  ## @param Jy_functor The functor depending on area coordinates for
  ## calculating the Jacobian determinant for \f$K_y\f$.
  ## @param sphere_center Center coordinates of the sphere manifold,
  ## which is stored in a row vector.
  ## @param sphere_radius The radius of the sphere model.
  ## @param ret The quadrature value.
  
  ## Generate the pullback of the kernel function \f$k(x, y)\f$ to the reference
  ## cell \f$\hat{K}_x \times \hat{K}_y\f$.
  k_loc = @(kx_area_coord, ky_area_coord) BEMKernelPullbackOnS2(kx_area_coord, ky_area_coord, nx_functor, ny_functor, kernel_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list, sphere_center, sphere_radius);
  
  ## Generate the pullback of \f$H_loc\f$ to the reference cell
  ## \f$\hat{K}_x \times \hat{K}_y\f$.
  H_loc = @(kx_area_coord, ky_area_coord) Jx_functor(kx_area_coord) * Jy_functor(ky_area_coord) * kx_basis_function(kx_area_coord) * ky_basis_function(ky_area_coord);
  
  ## Generate the integrand which combines both \f$H_{loc}\f$ and
  ## \f$K_{loc}\f$ on the reference cells.
  full_integrand_on_Khat = @(kx_area_coord, ky_area_coord) H_loc(kx_area_coord, ky_area_coord) * k_loc(kx_area_coord, ky_area_coord);
  
  ## Note: the full integrand above has its domain in the 2D unit simplex. It should further be pulled back to the master element \f$K^0 = \{(\xi_1, \xi_2) \in R^2 \big\vert 0 \leq \xi_1 \leq 1, 0 \leq \xi_2 \leq \xi_1\}\f$.
  full_integrand_on_K0 = @(kx_master_coord, ky_master_coord) full_integrand_on_Khat([1 - kx_master_coord(:, 1), kx_master_coord(:, 2)], [1 - ky_master_coord(:, 1), ky_master_coord(:, 2)]);
  
  ## Pullback the full integrand to the quadrature space \f$p = (\omega_1, \omega_2, \omega_3, \eta)\f$ and the integral domain is split into six parts.
  full_integrand_in_quad_space = @(p) p(:, 1) .* (
				   full_integrand_on_K0([sum(p(:, 1:3), 2), sum(p(:, 1:2), 2)],
							[p(:, 1) .* (1 - p(:, 4)) + p(:, 2) + p(:, 3), p(:, 2)]) +
				   full_integrand_on_K0([p(:, 1) .* (1 - p(:, 4)) + p(:, 2) + p(:, 3), p(:, 1) .* (1 - p(:, 4)) + p(:, 2)],
							[sum(p(:, 1:3), 2), p(:, 2)]) +
				   full_integrand_on_K0([sum(p(:, 1:3), 2), p(:, 1) .* p(:, 4) + p(:, 2)],
							[sum(p(:, 2:3), 2), p(:, 2)]) +
				   full_integrand_on_K0([p(:, 1) .* (1 - p(:, 4)) + p(:, 2) + p(:, 3), p(:, 2)],
							[sum(p(:, 1:3), 2), sum(p(:, 1:2), 2)]) +
				   full_integrand_on_K0([sum(p(:, 1:3), 2), p(:, 2)],
							[p(:, 1) .* (1 - p(:, 4)) + p(:, 2) + p(:, 3), p(:, 1) .* (1 - p(:, 4)) + p(:, 2)]) +
				   full_integrand_on_K0([sum(p(:, 2:3), 2), p(:, 2)],
							[sum(p(:, 1:3), 2), p(:, 1) .* p(:, 4) + p(:, 2)]));
  
  ## Generate the functor for the inner integral on a unit simplex \f$T^3\f$
  ## using tensor product based Gauss-Legendre quadrature rule.
  simplex_corners = GenerateSimplexCorners(3);
  integral_on_simplex = @(eta) tetra_tproduct(@(omega1, omega2, omega3) full_integrand_in_quad_space([omega1, omega2, omega3, eta]), norder_for_omega, simplex_corners(:, 1), simplex_corners(:, 2), simplex_corners(:, 3));
  
  ## Generate the 1D line using Gauss-Legendre rule.
  [qpts, qwts] = GaussLegendreRule(norder_for_eta);
  ## Adjust the Gauss-Legendre rule from \f$[-1, 1]\f$ to \f$[0, 1]\f$.
  [qpts, qwts] = rule_adjust(-1, 1, 0, 1, norder_for_eta, qpts, qwts);
  
  ret = ApplyQuadRule(integral_on_simplex, qpts', qwts');
endfunction
