function ret = Plot2DRect(varargin)
  ## Plot2DRect - Plot a 2D rectangle.
  ## @param bounding_box Bounding box of the rectangle, \f$[x_{\min}, y_{\min},
  ## w, h]\f$.
  ## @param line_style String specifying the line style.
  ## @param line_width Line width
  ## @param ret The figure handle.

  global fig_line_width;
  
  ret = [];

  if (length(varargin) >= 2)
    line_style = varargin{2};
  else
    line_style = "b-";
  endif

  if (length(varargin) >= 3)
    line_width = varargin{3};
  else
    line_width = fig_line_width;
  endif

  bounding_box = varargin{1};
  corner_list = [bounding_box(1), bounding_box(2); bounding_box(1) + bounding_box(3), bounding_box(2); bounding_box(1) + bounding_box(3), bounding_box(2) + bounding_box(4); bounding_box(1), bounding_box(2) + bounding_box(4); bounding_box(1), bounding_box(2)];

  ret = plot(corner_list(:, 1), corner_list(:, 2), line_style, "linewidth", line_width);
endfunction
