function PrintGCF(filename)
  ## PrintGCF - Print the current figure into files: eps and png.
  ## filename: base name of the figure file.

  h = gcf;
  print(h, cstrcat(filename, ".eps"), "-depsc2");
  ## print(h, cstrcat(filename, ".png"), "-dpng");
  system(["eps2png.sh " [filename, ".eps"]]);
endfunction
