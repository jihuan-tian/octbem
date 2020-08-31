function ret = JacobiP(x, n, alpha, beta)
  ## JacobiP - Jacobi polynomial.
  ## @param x A set of 1D coordinates.
  ## @param n Jacobi polynomial order.
  ## @param alpha \f$\alpha\f$ coefficient.
  ## @param beta \f$\beta\f$ coefficient.
  ## @param ret Polynomial values.

  switch (n)
    case 0
      ret = sqrt(2^(-alpha-beta-1) * gamma(alpha + beta + 2) / gamma(alpha + 1) / gamma(beta + 1)) * ones(size(x));
    case 1
      ret = 0.5 * JacobiP(x, 0, alpha, beta) * sqrt((alpha + beta + 3) / (alpha + 1) / (beta + 1)) .* ((alpha + beta + 2) * x + (alpha - beta));
    otherwise
      a_n = CalcA_n(n, alpha, beta);
      a_n_1 = CalcA_n(n - 1, alpha, beta);
      b_n_1 = CalcB_n(n - 1, alpha, beta);
      
      ret = (x .* JacobiP(x, n - 1, alpha, beta) - a_n_1 * JacobiP(x, n - 2, alpha, beta) - b_n_1 * JacobiP(x, n - 1, alpha, beta)) / a_n;
  endswitch
endfunction

function ret = CalcA_n(n, alpha, beta)
  ret = 2 / (2 * n + alpha + beta) * sqrt((n * (n + alpha + beta) * (n + alpha) * (n + beta)) / (2 * n + alpha + beta - 1) / (2 * n + alpha + beta + 1));
endfunction

function ret = CalcB_n(n, alpha, beta)
  ret = (beta^2 - alpha^2) / (2 * n + alpha + beta) / (2 * n + alpha + beta + 2);
endfunction
