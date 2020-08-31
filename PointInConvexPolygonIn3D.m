function point_to_polygon_relation = PointInConvexPolygonIn3D(p0, corner_nodes_of_polygon)
  ## PointInConvexPolygonIn3D - Determine if a point is inside, outside or on
  ## the boundary of a convex polygon in 3D space.
  ## @param p0 A list of target points.
  ## @param corner_nodes_of_polygon The list of ordered (clockwise or
  ## counter-clockwise) corner nodes for the polygon, which is stored as a N*3
  ## matrix.
  ## @param point_to_polygon_relation The relative relation between the point
  ## and the polygon: 1 for inside, 2 for outside, 3 for on the boundary.

  [number_of_targets, target_space_dim] = size(p0);
  [number_of_polygon_points, polygon_space_dim] = size(corner_nodes_of_polygon);

  if (number_of_polygon_points < 3)
    error("There should be more than three coplanar points defining the polygon!");
  endif

  if (polygon_space_dim != 3)
    error("The space dimension for the polygon should be 3!");
  endif

  if (target_space_dim != 3)
    error("The space dimension for the target points should be 3!");
  endif

  ## Calculate the projection distance of the point to the plane on which the
  ## polygon is defined.
  point_to_plane_distance = DistOfPointToPlaneIn3D(p0, corner_nodes_of_polygon);
  
  point_to_polygon_relation = zeros(number_of_targets, 1);
  for m = 1:number_of_targets
    if (abs(point_to_plane_distance(m)) < 1e-10)
      ## The point is already on the plane of the polygon.
      p0_to_edge_cross_products = zeros(number_of_polygon_points, 3);
      ## Iterate over each polygon edge to calculate the cross products.
      for n = 1:number_of_polygon_points
	if (n != number_of_polygon_points)
	  current_polygon_edge = corner_nodes_of_polygon(n + 1, :) - corner_nodes_of_polygon(n, :);
	else
	  current_polygon_edge = corner_nodes_of_polygon(1, :) - corner_nodes_of_polygon(n, :);
	endif
	p0_to_edge_cross_products(n, :) = NormalizeVect(cross(p0(m, :) - corner_nodes_of_polygon(n, :), current_polygon_edge));
      endfor

      ## Calculate the signed norms for all the cross product vectors in order
      ## to determine if the list of cross product vectors contains zero vector,
      ## which indicates the point lies on the polygon boundary. The direction
      ## of the first cross product vector is assumed to have the positive norm
      ## and used as reference.
      cross_product_signed_norms = zeros(number_of_polygon_points, 1);
      for k = 1:number_of_polygon_points
	cross_product_signed_norms(k) = sign(inner_prod(p0_to_edge_cross_products(k, :), p0_to_edge_cross_products(1, :))) * norm(p0_to_edge_cross_products(k, :));
      endfor
      
      if (length(find(abs(cross_product_signed_norms) < 1e-10)) > 0)
	## The point lies on the boundary of the polygon.
	point_to_polygon_relation(m) = 3;
      else
	## If the point does not lie on the boundary of the polygon, determine
	## if it lies inside or outside the polygon.
	if ((length(find(cross_product_signed_norms > 1e-10)) > 0)
	    && (length(find(cross_product_signed_norms < -1e-10)) > 0))
	  ## There exists both aligned and not aligned cross product vectors,
	  ## then the target point is outside the polygon.
	  point_to_polygon_relation(m) = 2;
	else
	  ## All cross product vectors are aligned in the same direction, then
	  ## the target point is inside the polygon.
	  point_to_polygon_relation(m) = 1;
	endif
      endif
    else
      point_to_polygon_relation(m) = 2;
    endif
  endfor
endfunction
