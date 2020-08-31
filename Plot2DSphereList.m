function ret = Plot2DSphereList(varargin)
  ## \ingroup OccProcSimu
  ## Plot2DSphereList - Plot a list of sphere objects.
  ## @param sphere_list A list of sphere objects stored in cell array.
  ## @param line_style String specifying the line style.
  ## @param line_width Line width.
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

  sphere_list = varargin{1};

  holdon_state = ishold();
  if (holdon_state == 0)
    hold on;
  endif

  for m = 1:length(sphere_list)
    ret = Plot2DCircle(sphere_list{m}.center, sphere_list{m}.radius, line_style, line_width);
  endfor

  if (holdon_state == 0)
    hold off;
  endif
endfunction
