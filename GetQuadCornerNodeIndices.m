function corner_node_indices = GetQuadCornerNodeIndices(order)
  ## GetQuadCornerNodeIndices - Get the indices for the corner nodes which are
  ## stored in the array of cell nodes with the specified order. Note: the
  ## ordering the quadrangle cell nodes depends on the implementation of the
  ## function AreaCoordsOnQuadrangle.
  ## @param order The specified order for representing the quadrangle's geometry.
  ## @param corner_node_indices The indices for the corner nodes.

  N = NumberOfNodesOnQuadrangle(order);
  corner_node_indices = [1, order + 1, N - order, N];
endfunction
