## This script print the linear Lagrange basis function using different color maps.
ConfigGraphicsToolkit;
pkg load image;

octave_color_map_names = {"viridis", "jet", "cubehelix", "hsv", "rainbow", "hot", "cool", "spring", "summer", "autumn", "winter", "gray", "bone", "copper", "pink", "ocean", "colorcube", "flag", "lines", "prism", "white"};

x_grid = linspace(-1, 1, 100);
y_grid = linspace(-1, 1, 100);
[xx, yy] = meshgrid(x_grid, y_grid);
zz = (xx - 1) .* (yy - 1) / 4;

## Plot into independent figures.
for m = 1:length(octave_color_map_names)
  ## Plot the original color map.
  figure(1);
  clf(1);
  
  colormap(octave_color_map_names{m});
  surf(xx, yy, zz, "EdgeColor", "none");
  view(2);
  grid off;
  colorbar;
  title(octave_color_map_names{m});
  
  PrintGCF(octave_color_map_names{m});

  ## Plot the reversed color map.
  figure(2);
  clf(2);
  
  colormap(ReverseColormap(octave_color_map_names{m}));
  surf(xx, yy, zz, "EdgeColor", "none");
  view(2);
  grid off;
  colorbar;
  title(cstrcat(octave_color_map_names{m}, "_reversed"));

  PrintGCF(cstrcat(octave_color_map_names{m}, "_reversed"));
endfor

## Plot into a single figure.
col_num = 4;
row_num = ceil((length(octave_color_map_names) - 1) * 2 / col_num);

figure(1);
clf(1);
for m = 1:(length(octave_color_map_names) - 1)
  ## Plot the original color map.
  subplot(row_num, col_num, m * 2 - 1);
  colormap(octave_color_map_names{m});
  subimage(ind2rgb(gray2ind(zz), colormap(octave_color_map_names{m})));
  grid off;
  axis off;
  axis equal;
  title(octave_color_map_names{m});

  ## Plot the reversed color map.
  subplot(row_num, col_num, m * 2);
  subimage(ind2rgb(gray2ind(zz), ReverseColormap(colormap(octave_color_map_names{m}))));
  grid off;
  axis off;
  axis equal;
  title(cstrcat(octave_color_map_names{m}, " reversed"));
endfor

print(gcf, cstrcat("color-maps", ".eps"), "-depsc2", "-S800,2400");
system(["eps2png.sh " ["color-maps", ".eps"]]);

close all;
