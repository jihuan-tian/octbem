function [qpts, qwts] = GaussJacobiLegendreRuleOnSquare(alpha, beta, orders)
  ## GaussLegendreJacobiRuleOnSquare - Generate the quadrature abscissas and
  ## weights in the 2D reference square: \f$[-1,1] \times [-1,1]\f$, where
  ## Gauss-Jacobi quadrature is used for the first dimension while
  ## Gauss-Legendre quadrature is used for the second dimension.
  ## @param alpha Exponent in \f$(1-x)^{\alpha}\f$ multiplied to the integrand
  ## for the second dimension.
  ## @param beta Exponent in \f$(1+x)^{\beta}\f$ multiplied to the integrand for
  ## the second dimension.
  ## @param orders The quadrature orders for the two dimensions.
  ## @param qpts Returned quadrature points.
  ## @param qwts Returned quadrature weights.

  if (length(orders) == 1)
    orders = orders * ones(1, 2);
  endif
  
  [qpts, qwts] = QuadRuleOnSquare(@(order) GaussJacobiRule(alpha, beta, order), @GaussLegendreRule, orders(1), orders(2));
endfunction
