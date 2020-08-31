function Plot3DAffineVectors(root_nodes, vector_nodes)
  ## Plot3DAffineVectors - Plot a list of 3D affine vectors.
  ## @param root_nodes A list of 3D root points which are stored in a Nx3 matrix.
  ## @param vector_nodes A list of 3D vectors which are stored in a Nx3 matrix.

  M = size(root_nodes, 1);
  N = size(vector_nodes, 1);

  if (M != N)
    error("The numbers of root nodes and vectors should be the same!");
  else
    if (N > 0)
      holdon_state = ishold();
      if (holdon_state == 0)
	hold on;
      endif

      for m = 1:M
	plot3([root_nodes(m, 1); root_nodes(m, 1) + vector_nodes(m, 1)], [root_nodes(m, 2); root_nodes(m, 2) + vector_nodes(m, 2)], [root_nodes(m, 3); root_nodes(m, 3) + vector_nodes(m, 3)], "marker", "o", "markerfacecolor", "r", "color", "g");
      endfor

      axis equal;

      if (holdon_state == 0)
	hold off;
      endif
    else
      error("There are no affine vectors to be plotted!");
    endif
  endif
endfunction
