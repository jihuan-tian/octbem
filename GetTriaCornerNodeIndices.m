function corner_node_indices = GetTriaCornerNodeIndices(order)
  ## GetTriaCornerNodeIndices - Get the indices for the corner nodes which are
  ## stored in the array of cell nodes with the specified order. Note: the
  ## ordering the triangle cell nodes depends on the implementation of the
  ## function AreaCoordsOnTriangle.
  ## @param order The specified order for representing the triangle's geometry.
  ## @param corner_node_indices The indices for the corner nodes.

  corner_node_indices = [1, order + 1, NumberOfNodesOnTriangle(order)];
endfunction
