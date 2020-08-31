function color_map_reversed = ReverseColormap(color_map)
  ## ReverseColormap - Reverse a color map.
  ## @param color_map The original color map.
  ## @param color_map_reversed The reversed color map.

  if ischar(color_map)
    color_map_reversed = flipud(eval(color_map));
  else
    color_map_reversed = flipud(color_map);
  endif
endfunction
