function palette = ColormapToMaxima(color_map_name, is_reversed)
  ## ColormapToMaxima - Export the color map to Maxima palette.
  ## @param color_map_name Name of the color map used in GNU Octave.
  ## @param is_reversed Should the map value be reversed, i.e. flip the color
  ## map matrix upside down (optional).
  ## @param palette String for palette represented in Maxima.

  if (!exist("is_reversed", "var"))
    is_reversed = false;
  endif
  
  if (is_reversed)
    color_map = ReverseColormap(color_map_name);
    color_map_name_suffix = "_reversed";
  else
    color_map = eval(color_map_name);
    color_map_name_suffix = "";
  endif

  palette = cstrcat("color_map_", color_map_name, color_map_name_suffix, " : [\"", strjoin(rgb2hex(color_map), "\",\""), "\"];");
endfunction
