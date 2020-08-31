function hex_str = rgb2hex(colors)
  ## rgb2hex - Converg a list of rgb colors to hex strings.

  color_num = size(colors, 1);
  hex_str = cell(color_num, 1);
  for m = 1:color_num
    hex_str{m} = cstrcat("#", sprintf("%02x", floor(colors(m, 1) * 255)), sprintf("%02x", floor(colors(m, 2) * 255)), sprintf("%02x", floor(colors(m, 3) * 255)));
  endfor
endfunction
