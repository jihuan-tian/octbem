function area_coord_matrix = AreaCoordsOnQuad(order)
  ## AreaCoordsOnQuad - Generate area coordinates on the reference
  ## quadrature with a specified order.
  ## @param order The specified order.
  ## @param area_coord_matrix The generated area coordinates in N*2 matrix.

  if (order >= 0)
    if (order == 0)
      area_coord_matrix = [0.5, 0.5];
    else
      subdiv_length = 1 / order;
      area_coord_matrix = zeros(NumberOfNodesOnQuadrangle(order), 2);
      counter = 1;

      ## Iterate over the first coordinate: \f$0 \leq m \leq 1\f$.
      for m = (0:subdiv_length:1)
	## Iterate over the second coordinate: \f$0 \leq n \leq 1 -
	## m\f$.
	for n = (0:subdiv_length:1)
	  area_coord_matrix(counter, :) = [m, n];
	  counter = counter + 1;
	endfor
      endfor
    endif
  else
    error("The order should be larger than or equal to 0.");
  endif
endfunction
