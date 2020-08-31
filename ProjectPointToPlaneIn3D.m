function p0_projection = ProjectPointToPlaneIn3D(p0, points_on_plane)
  ## ProjectPointToPlaneIn3D - Project a point to a plane in 3D space defined by
  ## a list of coplanar points.
  ## @param p0 A list of target points, which will be projected.
  ## @param points_on_plane The list of points defining the plane, which is
  ## stored as a N*3 matrix.
  ## @param p0_projection The projection of p0 onto the plane.

  [number_of_plane_points, plane_space_dim] = size(points_on_plane);
  [number_of_targets, target_space_dim] = size(p0);

  if (number_of_plane_points < 3)
    error("There should be more than three coplanar points defining the plane!");
  endif

  if (plane_space_dim != 3)
    error("The space dimension for the plane should be 3!");
  endif

  if (target_space_dim != 3)
    error("The space dimension for the target points should be 3!");
  endif
  
  [point_to_plane_distance, plane_normal_vector] = DistOfPointToPlaneIn3D(p0, points_on_plane);
  p0_projection = zeros(number_of_targets, 3);
  for m = 1:number_of_targets
    p0_projection(m, :) = p0(m, :) - plane_normal_vector * point_to_plane_distance(m);
  endfor
endfunction
