function ret = Plot2DPointList(varargin)
  ## Plot2DPointList - Plot a list of 2D points.
  ## @param point_list Matrix with dimension \f$N \times 2\f$ storing the point
  ## list. Each row represents a 2D point and the columns represents \f$XY\f$
  ## components.
  ## @param point_style String specifying the line type and/or point type for
  ## the plotting.
  ## @param ret The figure handle.

  global fig_line_width fig_marker_size;
  
  ret = [];

  if (length(varargin) >= 2)
    pt_style = varargin{2};
  else
    pt_style = "ro";
  endif
  
  ret = plot(varargin{1}(:,1), varargin{1}(:,2), pt_style);
endfunction
