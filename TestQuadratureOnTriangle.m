## This test code verifies that the quadrature formula for a triangle
## has the form \f$I = A \sum_i f(x_i) w_i\f$, where \f$A\f$ is the
## area of the triangle.

format long;

f = @(x) exp(-x(:, 2).^2) .* cos(x(:, 1) .* x(:, 2));
A = 0.5;

for order = 5:9
  [xtab, ytab, weight] = triangle_unit_product_set(order);
  ApplyQuadRule(f, [xtab', ytab'], weight') * A
endfor

[xtab, ytab, weight] = triangle_unit_set(6, 6);
ApplyQuadRule(f, [xtab', ytab'], weight') * A
