function [qpts, qwts] = GaussJacobiRuleOnCube(alphas, betas, orders)
  ## GaussJacobiRuleOnCube - Generate the Gauss-Jacobi quadrature abscissas and
  ## weights in the 3D reference cube: \f$[-1,1]^3\f$.
  ## @param alphas Exponents in \f$(1-x)^{\alpha}\f$ multiplied to the integrand
  ## for three dimensions.
  ## @param betas Exponents in \f$(1+x)^{\beta}\f$ multiplied to the integrand
  ## for three dimensions.
  ## @param orders The quadrature orders for the three dimensions.
  ## @param qpts Returned quadrature points.
  ## @param qwts Returned quadrature weights.

  if (length(alphas) == 1)
    alphas = alphas * ones(1, 3);
  endif

  if (length(betas) == 1)
    betas = betas * ones(1, 3);
  endif
  
  if (length(orders) == 1)
    orders = orders * ones(1, 3);
  endif

  [qpts, qwts] = QuadRuleOnCube(@(order) GaussJacobiRule(alphas(1), betas(1), order), @(order) GaussJacobiRule(alphas(2), betas(2), order), @(order) GaussJacobiRule(alphas(3), betas(3), order), orders(1), orders(2), orders(3));
endfunction
