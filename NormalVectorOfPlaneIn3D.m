function plane_normal = NormalVectorOfPlaneIn3D(coplanar_plane_points)
  ## NormalVectorOfPlaneIn3D - Calculate the normal vector of the plane defined by
  ## three coplanar points. The counter-clockwisely ordered nodes defines the
  ## positive direction of the normal vector.
  ## @param coplanar_plane_points Three coplanar points defining the plane.
  ## @param plane_normal The unit normal vector of the plane.

  [number_of_points, space_dim] = size(coplanar_plane_points);

  if (number_of_points < 3)
    error("There should be more than three coplanar points defining the plane!");
  endif

  if (space_dim != 3)
    error("The space dimension for the plane should be 3!");
  endif
  
  plane_normal = NormalizeVect(cross((coplanar_plane_points(2, :) - coplanar_plane_points(1, :)),
				     (coplanar_plane_points(3, :) - coplanar_plane_points(1, :))));
endfunction
