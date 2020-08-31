tria = [-2.499441205475897, 6.615921400291889, 2;
	-6.382273440256203, 3.662787464117798, -2;
	2.445576058240803, -6.253368887081748, 2.375194311054145];

## The point should be outside the triangle: 2
PointInConvexPolygonIn3D([-9.660004857787385, 7.212415070793341, 0], tria)
## The point should be outside the triangle: 2
PointInConvexPolygonIn3D([4.402836868910871, 8.503593533090886, 0], tria)
## The point should be inside the triangle: 1
PointInConvexPolygonIn3D([-2.763799942032662, 2.797088654296158, 0.701077676535478], tria)
## The point should be on the boundary of the triangle: 3
PointInConvexPolygonIn3D([-1.11841579641137, -2.250004254375943, 0.608834690880248], tria)
