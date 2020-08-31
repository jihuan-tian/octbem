function ret = ErichsenQuadCommonEdgeFlat(kernel_function, norder_for_eta,
					  norder_for_omega, kx_basis_function,
					  ky_basis_function,
					  kx_shape_functions_for_geometry,
					  ky_shape_functions_for_geometry,
					  kx_cell_node_coord_list,
					  ky_cell_node_coord_list,
					  nx, ny,
					  Jx, Jy)
  ## ErichsenQuadCommonEdge - Perform the common edge integration using the method in
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
  ## @param nx The cell normal vector for \f$K_x\f$.
  ## @param ny The cell normal vector for \f$K_y\f$.
  ## @param Jx The surface metric (Jacobian) for \f$K_x\f$ from reference to real cell.
  ## @param Jy The surface metric (Jacobian) for \f$K_y\f$ from reference to real cell.
  ## @param ret The quadrature value.
  
  ## Generate the pullback of the kernel function \f$k(x, y)\f$ to the reference
  ## cell \f$\hat{K}_x \times \hat{K}_y\f$.
  k_loc = @(kx_area_coord, ky_area_coord) BEMKernelPullbackFlat(kx_area_coord, ky_area_coord, nx, ny, kernel_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list);
  
  ## Generate the pullback of \f$H_loc\f$ to the reference cell
  ## \f$\hat{K}_x \times \hat{K}_y\f$.
  H_loc = @(kx_area_coord, ky_area_coord) Jx * Jy * kx_basis_function(kx_area_coord) * ky_basis_function(ky_area_coord);
  
  ## Generate the integrand which combines both \f$H_{loc}\f$ and \f$K_{loc}\f$ on the reference cells.
  full_integrand_on_Khat = @(kx_area_coord, ky_area_coord) H_loc(kx_area_coord, ky_area_coord) * k_loc(kx_area_coord, ky_area_coord);
  
  ## Note: the full integrand above has its domain in the 2D unit simplex. It should further be pulled back to the master element \f$K^0 = \{(\xi_1, \xi_2) \in \mathbb{R}^2 \big\vert 0 \leq \xi_1 \leq 1, 0 \leq \xi_2 \leq \xi_1\}\f$.
  full_integrand_on_K0 = @(kx_master_coord, ky_master_coord) full_integrand_on_Khat([1 - kx_master_coord(:, 1), kx_master_coord(:, 2)], [1 - ky_master_coord(:, 1), -ky_master_coord(:, 2)]);
  
  ## Pullback the full integrand to the Erichsen quadrature space \f$p = (\omega_1, \omega_2, \eta_1, \eta_2)\f$ and the integral domain is split into six parts.
  full_integrand_in_quad_space = @(p) p(:, 1).^2 .* p(:, 3) .* (
				   full_integrand_on_K0([sum(p(:, 1:2), 2), p(:, 1) .* p(:, 3) .* p(:, 4)],
							[p(:, 1) .* (1 - p(:, 3)) + p(:, 2), p(:, 1) .* (p(:, 3) - 1)]) +
				   full_integrand_on_K0([sum(p(:, 1:2), 2), p(:, 1) .* p(:, 3)],
							[p(:, 1) .* (1 - p(:, 4) .* p(:, 3)) + p(:, 2), p(:, 1) .* (p(:, 4) .* p(:, 3) - 1)]) +
				   full_integrand_on_K0([sum(p(:, 1:2), 2), p(:, 1)],
							[p(:, 1) .* p(:, 3) + p(:, 2), -p(:, 1) .* p(:, 4) .* p(:, 3)]) +
				   full_integrand_on_K0([p(:, 1) .* (1 - p(:, 3) + p(:, 4) .* p(:, 3)) + p(:, 2), p(:, 1) .* p(:, 4) .* p(:, 3)],
							[sum(p(:, 1:2), 2), -p(:, 1)]) +
				   full_integrand_on_K0([p(:, 1) .* p(:, 3) .* p(:, 4) + p(:, 2), p(:, 1) .* p(:, 3) .* p(:, 4)],
							[p(:, 1), -p(:, 1) .* p(:, 3)]) +
				   full_integrand_on_K0([p(:, 1) .* p(:, 3) + p(:, 2), p(:, 1) .* p(:, 3)],
							[sum(p(:, 1:2), 2), -p(:, 1) .* p(:, 4) .* p(:, 3)]));
  
  ## Generate Gauss-Legendre quadrature rule on the 2D simplex wrt omega.
  [simplex_qpts_omega1, simplex_qpts_omega2, simplex_qwts] = triangle_unit_product_set(norder_for_omega);
  
  ## Generate the functor for the inner integral on the 2D simplex \f$T^2\f$
  ## using tensor product based Gauss-Legendre quadrature rule.
  integral_on_simplex = @(etas) ApplyQuadRule(@(omegas) full_integrand_in_quad_space([omegas(1), omegas(2), etas(1), etas(2)]), [simplex_qpts_omega1', simplex_qpts_omega2'], simplex_qwts') * triangle_unit_volume();
  
  ## Integration on the 2D rectangle using Gauss-Legendre rule.
  [qpts_1d, qwts_1d] = GaussLegendreRule(norder_for_eta);
  ## Adjust the Gauss-Legendre rule from [-1, 1] to [0, 1].
  [qpts_1d, qwts_1d] = rule_adjust(-1, 1, 0, 1, norder_for_eta, qpts_1d, qwts_1d);
  qpts = tensor_prod2(qpts_1d, qpts_1d);
  qwts = tensor_prod2(qwts_1d, qwts_1d);
  
  ret = ApplyQuadRule(integral_on_simplex, qpts, qwts);
endfunction
