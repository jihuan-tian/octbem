function ret = LegendreP(x, n)
  ## LegendreP - Legendre polynomial.
  ## @param x A set of 1D coordinates.
  ## @param n Legendre polynomial order.
  ## @param ret Polynomial values.

  ret = JacobiP(x, n, 0, 0);
endfunction
