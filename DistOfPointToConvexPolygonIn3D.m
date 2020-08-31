function distance = DistOfPointToConvexPolygonIn3D(p0, corner_nodes_of_polygon)
  ## DistOfPointToConvexPolygonIn3D - Calculate the distance between a point
  ## and a convex polygon in 3D space.
  ## @param p0 The target point.
  ## @param corner_nodes_of_polygon The list of ordered (clockwise or
  ## counter-clockwise) corner nodes for the polygon, which is stored as a N*3
  ## matrix.
  ## @param The distance between the point and the polygon.

  [number_of_points, point_space_dim] = size(p0);
  [M, polygon_space_dim] = size(corner_nodes_of_polygon);

  if (point_space_dim != 3)
    error("The space dimension should be 3!");
  endif

  if (point_space_dim != polygon_space_dim)
    error("The space dimension for the point and the polygon should be the same!")
  endif

  ## Calculate the projection of the target point onto the plane of the polygon.
  x0 = ProjectPointToPlaneIn3D(p0, corner_nodes_of_polygon);
  ## Determine the position relation between the projection point and the
  ## polygon.
  x0_to_polygon_relation = PointInConvexPolygonIn3D(x0, corner_nodes_of_polygon);

  distance = zeros(number_of_points, 1);
  for m = 1:number_of_points
    if (x0_to_polygon_relation(m) == 2)
      ## If the projection point x0 is outside the polygon, calculate the
      ## distance between each polygon node and the target point, then select
      ## the minimum as the point to polygon distance.
      point_to_polygon_node_distance = zeros(M, 1);
      for n = 1:M
	point_to_polygon_node_distance(n) = norm(p0(m, :) - corner_nodes_of_polygon(n, :));
      endfor
      distance(m) = min(point_to_polygon_node_distance);
    else
      ## If the projection point x0 is inside or on the boundary of the polygon,
      ## calculate the distance between the target point and the plane of the
      ## polygon.
      distance(m) = abs(DistOfPointToPlaneIn3D(p0(m, :), corner_nodes_of_polygon));
    endif
  endfor
endfunction
