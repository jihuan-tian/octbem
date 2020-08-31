function [qpts, qwts] = QuadRuleOnSquare(quad_rule_x, quad_rule_y, order_x, order_y)
  ## QuadRuleOnSquare - Generate the quadrature abscissas and weights in the 2D
  ## reference square: \f$[-1,1]^2\f$ with specified rules for two dimensions.
  ## @param quad_rule_x The quadrature rule for the first dimension.
  ## @param quad_rule_y The quadrature rule for the second dimension.
  ## @param order_x The quadrature order in the first dimension.
  ## @param order_y The quadrature order in the second dimension.
  ## @param qpts Returned quadrature points.
  ## @param qwts Returned quadrature weights.

  [x_qpts, x_qwts] = quad_rule_x(order_x);
  [y_qpts, y_qwts] = quad_rule_y(order_y);

  qpts = tensor_prod2(x_qpts, y_qpts);
  qwts = tensor_prod2(x_qwts, y_qwts);
endfunction
