function [qpts, qwts] = GaussLegendreRuleOnSquare(orders)
  ## GaussLegendreRuleOnSquare - Generate the Gauss-Legendre quadrature
  ## abscissas and weights in the 2D reference square: \f$[-1,1] \times
  ## [-1,1]\f$.
  ## @param orders The Gauss-Legendre quadrature orders for the two dimensions.
  ## @param qpts Returned quadrature points.
  ## @param qwts Returned quadrature weights.

  if (length(orders) == 1)
    orders = orders * ones(1, 2);
  endif
  
  [qpts, qwts] = QuadRuleOnSquare(@GaussLegendreRule, @GaussLegendreRule, orders(1), orders(2));
endfunction
