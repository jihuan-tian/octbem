function [qpts, qwts] = GaussLegendreRuleOnCube(orders)
  ## GaussLegendreRuleOnCube - Generate the Gauss-Legendre quadrature abscissas
  ## and weights in the 3D reference cube: \f$[-1,1]^3\f$.
  ## @param orders The quadrature orders for the three dimensions.
  ## @param qpts Returned quadrature points.
  ## @param qwts Returned quadrature weights.

  if (length(orders) == 1)
    orders = orders * ones(1, 3);
  endif

  [qpts, qwts] = QuadRuleOnCube(@GaussLegendreRule, @GaussLegendreRule, @GaussLegendreRule, orders(1), orders(2), orders(3));
endfunction
