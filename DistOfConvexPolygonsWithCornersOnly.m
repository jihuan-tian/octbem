function distance = DistOfConvexPolygonsWithCornersOnly(corner_nodes_of_polygon1,
							corner_nodes_of_polygon2)
  ## DistOfConvexPolygonsWithCornersOnly - Calculate the distance between two
  ## convex polygons which may not be coplanar in 3D space. Each polygon is
  ## represented by a list of corner node coordinates.
  ## @param corner_nodes_of_polygon1 The list of ordered (clockwise or
  ## counter-clockwise) corner nodes for the first polygon, which is stored as a
  ## N*3 matrix.
  ## @param corner_nodes_of_polygon1 The list of ordered (clockwise or
  ## counter-clockwise) corner nodes for the second polygon, which is stored as
  ## a N*3 matrix.
  ## @param distance The distance between the two polygons.

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

  ## Calculate the distance from nodes in polygon 1 to polygon 2.
  polygon1_nodes_to_polygon2_distance = DistOfPointToConvexPolygonIn3D(corner_nodes_of_polygon1, corner_nodes_of_polygon2);
  ## Calculate the distance from the center of polygon 1 to polygon 2.
  polygon1_center = sum(corner_nodes_of_polygon1) / M;
  polygon1_center_to_polygon2_distance = DistOfPointToConvexPolygonIn3D(polygon1_center, corner_nodes_of_polygon2);
  
  ## Calculate the distance from nodes in polygon 2 to polygon 1.
  polygon2_nodes_to_polygon1_distance = DistOfPointToConvexPolygonIn3D(corner_nodes_of_polygon2, corner_nodes_of_polygon1);
  ## Calculate the distance from the center of polygon 2 to polygon 1.
  polygon2_center = sum(corner_nodes_of_polygon2) / N;
  polygon2_center_to_polygon1_distance = DistOfPointToConvexPolygonIn3D(polygon2_center, corner_nodes_of_polygon1);

  distance = min([polygon1_nodes_to_polygon2_distance; polygon1_center_to_polygon2_distance; polygon2_nodes_to_polygon1_distance; polygon2_center_to_polygon1_distance]);
endfunction
