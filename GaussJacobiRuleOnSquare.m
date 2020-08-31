function [qpts, qwts] = GaussJacobiRuleOnSquare(alphas, betas, orders)
  ## GaussJacobiRuleOnSquare - Generate the Gauss-Jacobi quadrature abscissas
  ## and weights in the 2D reference square: \f$[-1,1] \times [-1,1]\f$.
  ## @param alphas Exponents in \f$(1-x)^{\alpha}\f$ multiplied to the integrand
  ## for two dimensions.
  ## @param betas Exponents in \f$(1+x)^{\beta}\f$ multiplied to the integrand
  ## for two dimensions.
  ## @param orders The Gauss-Jacobi quadrature orders for the two dimensions.
  ## @param qpts Returned quadrature points.
  ## @param qwts Returned quadrature weights.

  if (length(alphas) == 1)
    alphas = alphas * ones(1, 2);
  endif

  if (length(betas) == 1)
    betas = betas * ones(1, 2);
  endif

  if (length(orders) == 1)
    orders = orders * ones(1, 2);
  endif
  
  [qpts, qwts] = QuadRuleOnSquare(@(order) GaussJacobiRule(alphas(1), betas(1), order), @(order) GaussJacobiRule(alphas(2), betas(2), order), orders(1), orders(2));
endfunction
