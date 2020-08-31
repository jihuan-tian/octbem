function [norder_for_eta, norder_for_omega] = ErichsenSingularQuadOrder(basis_function_polynomial_order, function_space_sobolev_order, mesh_size_estimate)
  ## ErichsenSingularQuadOrder - Calculate the Erichsen quadrature order
  ## (actually, the number of quadrature points) for singular integrals. The
  ## integration space is \f$\eta \times \omega \in R_d \times T_{4 - d}\f$,
  ## where \f$\omega\f$ is firstly integreated, then \f$\eta\f$.
  ## @param basis_function_polynomial_order The polynomial order \f$p\f$ of the
  ## adopted basis function.
  ## @param function_space_sobolev_order The Sobolev space order for the adopted
  ## function space, i.e. \f$s_1\f$
  ## @param mesh_size_estimate An estimate of the mesh size \f$h\f$.
  ## @param norder_for_eta Number of quadrature points for the integration
  ## variable \f$\eta\f$.
  ## @param norder_for_omega Number of quadrature points for the integration
  ## variable \f$\omega\f$.
  
  norder_for_eta = ceil((basis_function_polynomial_order - function_space_sobolev_order) / 2 * abs(log(mesh_size_estimate)));
  norder_for_omega = ceil(basis_function_polynomial_order - function_space_sobolev_order);
endfunction
