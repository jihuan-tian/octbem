function order = QuadGeometryOrder(number_of_nodes)
  ## QuadGeometryOrder - Calculate the geometry representation order of the
  ## reference quadrangle by deriving from the number of reference nodes. This
  ## function is the inverse of NumberOfNodesOnQuad.
  ## @param number_of_nodes Total number of nodes.
  ## @param order The geometry representation order.

  order = sqrt(number_of_nodes) - 1;
endfunction
