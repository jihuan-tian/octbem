function h = Cell3DRange(cell_corners)
  ## Cell3DRange - Calculate the range of the cell in 3D space, which is the
  ## maximum inter-distance between corner nodes.
  ## @param cell_corners The 3D coordinates for cell corners, stored as a n*3 matrix.
  ## @param h The cell range, which is the maximum inter-distance between corner nodes.

  N = size(cell_corners, 1);
  number_of_node_pairs = N * (N - 1) / 2;
  inter_dists = zeros(number_of_node_pairs, 1);
  
  count = 0;
  for m = 2:N
    for n = 1:(m - 1)
      count += 1;
      inter_dists(count) = norm(cell_corners(m, :) - cell_corners(n, :));
    endfor
  endfor

  h = max(inter_dists);
endfunction
