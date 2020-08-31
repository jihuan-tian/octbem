## This script convert all Octave colormaps to Maxima.
octave_color_map_names = {"viridis", "jet", "cubehelix", "hsv", "rainbow", "hot", "cool", "spring", "summer", "autumn", "winter", "gray", "bone", "copper", "pink", "ocean", "colorcube", "flag", "lines", "prism", "white"};

## Export the color maps.
strjoin(cellfun(@ColormapToMaxima, octave_color_map_names, "UniformOutput", false), "\n")
## Export the reversed color maps.
strjoin(cellfun(@(X) ColormapToMaxima(X, true), octave_color_map_names, "UniformOutput", false), "\n")
