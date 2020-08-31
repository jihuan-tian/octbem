function ret = Plot2DEllipse(varargin)
  ## \ingroup RadarSimu
  ## PlotEllipse - Plot a 2D ellipse.
  ## @param center_coord: \f$(x,y)\f$ center coordinate.
  ## @param long_axis_dir Long axis direction.
  ## @param a Half length of the longer axis.
  ## @param b Half length of the shorter axis.
  ## @param angle_range Angle range in radian.
  ## @param enable_draw_axis Enable drawing of axis.
  ## @param line_style String Specifying the line style.
  ## @param line_width Line width.
  ## @param ret The figure handle.

  global fig_line_width;
  
  ret = [];

  if (length(varargin) >= 5)
    angle_range = varargin{5};
  else
    angle_range = [0, 2 * pi];
  endif

  if (length(varargin) >= 6)
    enable_drawing_axis = varargin{6};
  else
    enable_drawing_axis = 0;
  endif
  
  if (length(varargin) >= 7)
    line_style = varargin{7};
  else
    line_style = "b-";
  endif

  if (length(varargin) >= 8)
    line_width = varargin{8};
  else
    line_width = fig_line_width;
  endif

  center_coord = varargin{1};
  center_coord = reshape(center_coord, 2, 1);
  long_axis_dir = varargin{2};
  long_axis_dir = reshape(long_axis_dir, 2, 1);
  long_axis_dir = long_axis_dir / norm(long_axis_dir);
  short_axis_dir = Rz2D(pi / 2) * long_axis_dir;

  a = varargin{3};
  b = varargin{4};

  dtheta = 0.001;
  ## Generate the standard ellipse centered at origin with no rotation.
  theta = linspace(angle_range(1), angle_range(2), floor((angle_range(2) - angle_range(1)) / dtheta));
  ellipse_points = [a * cos(theta); b * sin(theta)];
  ## Rotation angle of the ellipse.
  rot_angle = atan2(long_axis_dir(2), long_axis_dir(1));
  ## Rotate and translat the standard ellipse.
  ellipse_points = Rz2D(rot_angle) * ellipse_points + center_coord;
  
  ret = plot(ellipse_points(1, :), ellipse_points(2, :), line_style, "linewidth", line_width);

  if (enable_drawing_axis)
    holdon_state = ishold();
    if (holdon_state == 0)
      hold on;
    endif

    ## Plot the long axis.
    long_axis_end1 = center_coord + a * long_axis_dir;
    long_axis_end2 = center_coord - a * long_axis_dir;
    ret = line([long_axis_end1(1), long_axis_end2(1)], [long_axis_end1(2), long_axis_end2(2)], "linestyle", "--", "color", "b", "linewidth", 1);
    ## Plot the short axis.
    short_axis_end1 = center_coord + b * short_axis_dir;
    short_axis_end2 = center_coord - b * short_axis_dir;
    ret = line([short_axis_end1(1), short_axis_end2(1)], [short_axis_end1(2), short_axis_end2(2)], "linestyle", "--", "color", "b", "linewidth", 1);
    
    if (holdon_state == 0)
      hold off;
    endif
  endif
endfunction
