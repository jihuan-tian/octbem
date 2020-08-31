function number_of_nodes = NumberOfNodesOnTriangle(order)
  ## NumberOfNodesOnTriangle - Calculate the total number of nodes to be
  ## distributed on the reference triangle with the specified geometry
  ## representation order. This functions is the inverse of
  ## TriaGeometryOrder.
  ## @param order The specified order for representing the triangle's geometry.
  ## @param number_of_nodes Total number of nodes.

  number_of_nodes = (order + 2) * (order + 1) / 2;
endfunction
