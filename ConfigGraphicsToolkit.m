## Declare global variables used for graphics configuration.
DeclareGraphicsGlobals;

## Select default graphical toolkit.
graphics_toolkit("gnuplot");

## Set default plotting parameters.
fig_font_size = 14;
legend_font_size = 8;
fig_line_width = 3;
fig_marker_size = 6;
set (0, "defaultaxesfontname", "Helvetica");
set (0, "defaultaxesfontsize", fig_font_size);
set (0, "defaulttextfontname", "Helvetica");
set (0, "defaulttextfontsize", fig_font_size);

## Generate a list of colors used for line plotting.
clist = GenColorList();

## The initial figure handle.
h = 0;
