function i = function_2d_set ( action, i )

%*****************************************************************************80
%
%% FUNCTION_2D_SET sets or reports the index of the current 2D function.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    07 April 2008
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character ACTION(*):
%    'GET', please return the index of the current function.
%    'SET', please set the function according to the input index.
%
%    Input/output, integer I.
%    If ACTION = 'SET', then I is input, and is the index of the desired
%    function.
%    If ACTION = 'GET', then I is output, and is the index of the current
%    function.
%
  persistent i_save;

  if ( action == 'SET' )
    i_save = i;
  elseif ( action == 'GET' )
    i = i_save;
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, 'FUNCTION_2D_SET - Fatal error!\n' );
    fprintf ( 1, '  Unrecognized input argument ACTION = "%s"\n', action );
    error ( 'FUNCTION_2D_SET - Fatal error!' );
  end

  return
end
