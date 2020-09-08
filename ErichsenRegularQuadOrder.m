function norder = ErichsenRegularQuadOrder(kernel_singularity_order, basis_function_polynomial_order, function_space_sobolev_order, mesh_size_estimate, panel_distance, galerkin_estimate_norm_index)
  ## ErichsenRegularQuadOrder - Calculate the Erichsen quadrature order
  ## (actually, the number of quadrature points) for regular integrals.
  ## @param kernel_singularity_order Singularity order of the kernel to be
  ## integrated.
  ## @param basis_function_polynomial_order The polynomial order \f$p\f$ of the
  ## adopted basis function.
  ## @param function_space_sobolev_order The Sobolev space order for the adopted
  ## function space, i.e. \f$s_1\f$
  ## @param mesh_size_estimate An estimate of the mesh size \f$h\f$.
  ## @param panel_distance The distance \f$\delta\f$ between the two panels.
  ## @param galerkin_estimate_norm_index The index of the Sobolev norm, which is
  ## used for measuring the error in Galerkin esitmate. It is equal to \f$s_1 -
  ## t\f$. When it is zero, the error is \f$L^2\f$.
  ## @param norder Number of quadrature points.

  
  # \f$t\f$ should be range \f$[0, p+1-s_1]\f$.
  t = 0;
  D = max([0, log(panel_distance) / log(mesh_size_estimate / panel_distance)]);
  norder = ceil(0.5 * (2 * basis_function_polynomial_order + t - function_space_sobolev_order + (kernel_singularity_order + basis_function_polynomial_order + t - function_space_sobolev_order) * D));

  ## Limit the norder.
  max_norder = 5;
  if (norder > max_norder)
    norder = max_norder;
  endif
endfunction
