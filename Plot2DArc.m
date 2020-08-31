function ret = Plot2DArc(varargin)
  ## \ingroup RadarSimu
  ## Plot2DArc - Plot a 2D arc.
  ## @param center_coord \f$(x,y)\f$ center coordinate.
  ## @param r Radius.
  ## @param angle_range Angle range in radian.
  ## @param line_style String specifying the line style.
  ## @param line_width Line width.
  ## @param ret The figure handle.

  global fig_line_width;
  
  ret = [];

  if (length(varargin) >= 4)
    line_style = varargin{4};
  else
    line_style = "b-";
  endif

  if (length(varargin) >= 5)
    line_width = varargin{5};
  else
    line_width = fig_line_width;
  endif

  center_coord = varargin{1};
  r = varargin{2};
  angle_range = varargin{3};

  dtheta = 0.001;
  theta = linspace(angle_range(1), angle_range(2), floor((angle_range(2) - angle_range(1)) / dtheta));
  ret = plot(r * cos(theta) + center_coord(1), r * sin(theta) + center_coord(2), line_style, "linewidth", line_width);
endfunction
