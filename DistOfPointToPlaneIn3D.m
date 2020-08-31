function [distance, plane_normal_vector] = DistOfPointToPlaneIn3D(p0, points_on_plane)
  ## DistOfPointToPlaneIn3D - Calculate the signed distance between a point and
  ## a plane. They may not be coplanar.
  ## @param p0 A list of target points.
  ## @param points_on_plane The list of points defining the plane, which is
  ## stored as a N*3 matrix.
  ## @param distance The signed distance between the point and the plane.
  ## @param plane_normal_vector The unit normal vector of the plane.
  
  ## Reference: http://mathworld.wolfram.com/Point-PlaneDistance.html

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
  
  ## Use the first three points from the point list to define the plane, i.e.
  ## calculate the normal vector of the plane.
  plane_normal_vector = NormalVectorOfPlaneIn3D(points_on_plane(1:3, :));

  distance = zeros(number_of_targets, 1);
  for m = 1:number_of_targets
    ## Project the difference vector between the target point and any point in
    ## the list to the normal vector of the plane.
    distance(m) = inner_prod(p0(m, :) - points_on_plane(1, :), plane_normal_vector);
  endfor
endfunction
