function ret = ApplyQuadRule(functor, qpts, qwts, J)
  ## ApplyQuadRule - Apply a quadrature rule to a function.
  ## @param functor The function to be integrated, which should have been pulled
  ## back onto the reference cell in which the quadrature points are generated.
  ## The functor should accept a vector \f$x\f$ as its argument, instead of a
  ## list of coordinate components.
  ## @param qpts A collection of quadrature points, with dimension \f$N d\f$,
  ## where \f$N\f$ is the number of points, \f$d\f$ is the space dimension.
  ## @param qwts A collection of quadrature weights, , with dimension \f$N d\f$,
  ## where \f$N\f$ is the number of points, \f$d\f$ is the space dimension.
  ## @param J The Jacobian, which is a functor, maps from the reference cell to
  ## the real cell on which the functor originally is defined.

  ret = 0;
  N = size(qpts, 1);

  if (!exist("J", "var"))
    J = @(x) 1;
  endif
  
  for m = 1:N
    ## The product of the Jacobian and the quadrature weights.
    JxW = J(qpts(m, :)) * prod(qwts(m, :), 2);
    ret = ret + functor(qpts(m, :)) * JxW;
  endfor
endfunction
