function distance = DistOfCoplanarPointToConvexPolygonIn3D(p0, corner_nodes_of_polygon)
  ## DistOfCoplanarPointToConvexPolygonIn3D - Calculate the distance between a point
  ## and a coplanar convex polygon in 3D space.
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

  distance = zeros(number_of_points, 1);
  for k = 1:number_of_points
    point_to_polygon_distance = zeros(1, M);
    for m = 1:M			# Loop over each edge of the polygon.
      ## Get the 1st ending node of the current edge in the polygon.
      p1 = corner_nodes_of_polygon(m, :);
      ## Get the 2nd ending node of the current edge in the polygon.
      if (m != M)
	p2 = corner_nodes_of_polygon(m + 1, :);
      else
	p2 = corner_nodes_of_polygon(1, :);
      endif

      ## Direction vector of the edge.
      p1p2_dir = NormalizeVect(p2 - p1);
      ## Project \f$\overline{P_1 P_0}\f$ onto the direction of \f$\overline{P_1 P_2}\f$
      projection_of_vertical_line = inner_prod(p0(k, :) - p1, p1p2_dir);
      ## Intersection point of the vertical line and the edge.
      x0 = p1 + projection_of_vertical_line * p1p2_dir;

      if ((projection_of_vertical_line < 0) ||
	  (projection_of_vertical_line > norm(p2 - p1)))
	## The intersection point of the vertical line and the edge lies outside
	## the edge segment.
	point_to_polygon_distance(m) = min([norm(p0(k, :) - p1), norm(p0(k, :) - p2)]);
      else
	## The intersection point of the vertical line and the edge lies within
	## the edge segment.
	point_to_polygon_distance(m) = norm(p0(k, :) - x0);
      endif
    endfor

    distance(k) = min(point_to_polygon_distance);
  endfor
endfunction
