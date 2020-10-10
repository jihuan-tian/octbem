function number_of_nodes = NumberOfNodesOnQuadrangle(order)
  ## NumberOfNodesOnQuadrangle - Calculate the total number of nodes to be
  ## distributed on the reference quadrangle with the specified geometry
  ## representation order. This functions is the inverse of
  ## QuadGeometryOrder.
  ## @param order The specified order for representing the quadrangle's geometry.
  ## @param number_of_nodes Total number of nodes.

  number_of_nodes = (order + 1)^2;
endfunction
