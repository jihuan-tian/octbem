function distance = DistOfCoplanarConvexPolygons(corner_nodes_of_polygon1,
						 corner_nodes_of_polygon2)
  ## DistOfConvexPolygons - Calculate the distance between two coplanar convex
  ## polygons in 3D space. Each polygon is represented by a list of corner node
  ## coordinates.
  ## @param corner_nodes_of_polygon1 The list of ordered (clockwise or
  ## counter-clockwise) corner nodes for the first polygon, which is stored as a
  ## N*3 matrix.
  ## @param corner_nodes_of_polygon1 The list of ordered (clockwise or
  ## counter-clockwise) corner nodes for the second polygon, which is stored as
  ## a N*3 matrix.
  ## @param distance The distance between the two coplanar polygons.

  ## Number of corner nodes in polygon 1 and space dimension.
  [M, space_dim_1] = size(corner_nodes_of_polygon1);
  ## Number of corner nodes in polygon 2 and space dimension.
  [N, space_dim_2] = size(corner_nodes_of_polygon2);

  if (space_dim_1 != 3)
    error("The space dimension should be 3!");
  endif

  if (space_dim_1 != space_dim_2)
    error("The space dimension for the two polygons should be the same!")
  endif

  if (M < 3)
    error("There should be 3 points or more corners for the first polygon!");
  endif

  if (N < 3)
    error("There should be 3 points or more corners for the second polygon!");
  endif

  ## Distance from each node in polygon 1 to polygon 2.
  node_in_polygon1_to_polygon2_distance = zeros(1, M);
  
  for m = 1:M			# Loop over each node of polygon 1.
    node_in_polygon1_to_polygon2_distance(m) = DistOfCoplanarPointToConvexPolygonIn3D(corner_nodes_of_polygon1(m, :), corner_nodes_of_polygon2);
  endfor

  ## Distance from each node in polygon 2 to polygon 1.
  node_in_polygon2_to_polygon1_distance = zeros(1, N);
  
  for n = 1:N			# Loop over each node of polygon 2.
    node_in_polygon2_to_polygon1_distance(n) = DistOfCoplanarPointToConvexPolygonIn3D(corner_nodes_of_polygon2(n, :), corner_nodes_of_polygon1);
  endfor

  distance = min([min(node_in_polygon1_to_polygon2_distance), min(node_in_polygon2_to_polygon1_distance)]);
endfunction
