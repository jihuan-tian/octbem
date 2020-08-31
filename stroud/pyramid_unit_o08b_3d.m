function result = pyramid_unit_o08b_3d ( func )

%*****************************************************************************80
%
%% PYRAMID_UNIT_O08B_3D approximates an integral inside the unit pyramid in 3D.
%
%  Discussion:
%
%    An 8 point formula is used.
%
%    The (X,Y,Z) integration region can be represented as:
%
%    - ( 1 - Z ) <= X <= 1 - Z
%    - ( 1 - Z ) <= Y <= 1 - Z
%              0 <= Z <= 1.
%
%    When Z is zero, the integration region is a square lying in the (X,Y)
%    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
%    radius of the square diminishes, and when Z reaches 1, the square has
%    contracted to the single point (0,0,1).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    01 April 2008
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Carlos Felippa,
%    A compendium of FEM integration formulas for symbolic work,
%    Engineering Computation,
%    Volume 21, Number 8, 2004, pages 867-890.
%
%  Parameters:
%
%    Input, external FUNC, the name of the user supplied function which
%    evaluates F(X,Y,Z), of the form
%      function value = func ( x, y, z )
%
%    Output, real RESULT, the approximate integral of the function.
%
  order = 8;

  w = [ ...
   0.16438287736328777572, ...
   0.16438287736328777572, ...
   0.16438287736328777572, ...
   0.16438287736328777572, ...
   0.085617122636712224276, ...
   0.085617122636712224276, ...
   0.085617122636712224276, ...
   0.085617122636712224276 ]';
  x = [ ...
  -0.51197009372656270107, ...
   0.51197009372656270107, ...
   0.51197009372656270107, ...
  -0.51197009372656270107, ...
  -0.28415447557052037456, ...
   0.28415447557052037456, ...
   0.28415447557052037456, ...
  -0.28415447557052037456 ]';
  y = [ ...
  -0.51197009372656270107, ...
  -0.51197009372656270107, ...
   0.51197009372656270107, ...
   0.51197009372656270107, ...
  -0.28415447557052037456, ...
  -0.28415447557052037456, ...
   0.28415447557052037456, ...
   0.28415447557052037456 ]';
  z = [ ...
   0.11024490204163285720, ...
   0.11024490204163285720, ...
   0.11024490204163285720, ...
   0.11024490204163285720, ...
   0.518326526529795714229, ...
   0.518326526529795714229, ...
   0.518326526529795714229, ...
   0.518326526529795714229 ]';
%
%  Quadrature.
%
  quad = 0.0;
  for i = 1 : order
    quad = quad + w(i) * feval ( func, x(i), y(i), z(i) );
  end
%
%  Volume.
%
  volume = pyramid_unit_volume_3d ( );
%
%  Result.
%
  result = quad * volume;

  return
end
