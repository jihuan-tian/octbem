function [qpts, qwts] = QuadRuleOnCube(quad_rule_x, quad_rule_y, quad_rule_z, order_x, order_y, order_z)
  ## QuadRuleOnCube - Generate the quadrature abscissas and weights in the 3D
  ## reference cube: \f$[-1,1]^3\f$ with specified rules for three dimensions.
  ## @param quad_rule_x The quadrature rule for the first dimension.
  ## @param quad_rule_y The quadrature rule for the second dimension.
  ## @param quad_rule_z The quadrature rule for the third dimension.
  ## @param order_x The quadrature order in the first dimension.
  ## @param order_y The quadrature order in the second dimension.
  ## @param order_z The quadrature order in the third dimension.
  ## @param qpts Returned quadrature points.
  ## @param qwts Returned quadrature weights.

  ## Generate the quadrature points and weights for each dimension.
  [x_qpts, x_qwts] = quad_rule_x(order_x);
  [y_qpts, y_qwts] = quad_rule_y(order_y);
  [z_qpts, z_qwts] = quad_rule_z(order_z);

  qpts = tensor_prod3(x_qpts, y_qpts, z_qpts);
  qwts = tensor_prod3(x_qwts, y_qwts, z_qwts);
endfunction
