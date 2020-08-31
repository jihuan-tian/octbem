function [ o, x, w ] = epn_glg_01_1 ( n, alpha )

%*****************************************************************************80
%
%% EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
%
%  Discussion:
%
%    The rule has order O = 1.
%
%    The rule has precision P = 1.
%
%    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
%    Laguerre weight function:
%
%      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    29 January 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the spatial dimension.
%
%    Input, real ALPHA, the exponent of X in the weight function.
%    -1.0 < ALPHA.
%
%    Input, integer O, the order.
%
%    Output, real X(N,O), the abscissas.
%
%    Output, real W(O), the weights.
%
  if ( alpha <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'EPN_GLG_01_1 - Fatal error!\n' );
    fprintf ( 1, '  ALPHA <= -1.0\n' );
    error ( 'EPN_GLG_01_1 - Fatal error!' );
  end

  expon = 0;
  value1 = ep1_glg_monomial_integral ( expon, alpha );
  volume = value1 ^ n;

  expon = 1;
  value2 = ep1_glg_monomial_integral ( expon, alpha );

  o = 1;

  x = zeros ( n, o );
  w = zeros ( o, 1 );

  k = 0;
%
%  1 point.
%
  k = k + 1;
  x(1:n,k) = value2 / value1;
  w(k) = volume;

  return
end
