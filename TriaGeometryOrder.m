function order = TriaGeometryOrder(number_of_nodes)
  ## TriaGeometryOrder - Calculate the geometry representation order of the
  ## reference triangle by deriving from the number of reference nodes. This
  ## function is the inverse of NumberOfNodesOnTriangle.
  ## @param number_of_nodes Total number of nodes.
  ## @param order The geometry representation order.

  order = (sqrt(8 * number_of_nodes + 1) - 3) / 2;
endfunction
