function ret = LegendrePN(x, n)
  ## LegendrePN - Normalized Legendre polynomial.
  ## @param x A set of 1D coordinates.
  ## @param n Legendre polynomial order.
  ## @param ret Polynomial values.

  ret = LegendreP(x, n) / (0.5 + n)^0.5;
endfunction
