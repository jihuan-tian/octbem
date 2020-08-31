function ret = Plot2DCircle(varargin)
  ## \ingroup RadarSimu
  ## Plot2DCircle - Plot a 2D circle.
  ## @param center_coord \f$(x,y)\f$ center coordinate.
  ## @param r Radius.
  ## @param line_style String specifying the line style.
  ## @param line_width Line width.
  ## @param ret The figure handle.

  global fig_line_width;
  
  ret = [];

  if (length(varargin) >= 3)
    line_style = varargin{3};
  else
    line_style = "b-";
  endif

  if (length(varargin) >= 4)
    line_width = varargin{4};
  else
    line_width = fig_line_width;
  endif

  center_coord = varargin{1};
  r = varargin{2};

  theta = linspace(0, 2*pi, 3600);
  ret = plot(r * cos(theta) + center_coord(1), r * sin(theta) + center_coord(2), line_style, "linewidth", line_width);
endfunction
