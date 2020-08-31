## Original 4th order node coordinates for a flat triangle.
tria_node_global_coord_flat = [0, 0, 5;
			       0.43, 0.58, 3.75;
			       0.86, 1.15, 2.5;
			       1.29, 1.73, 1.25;
			       1.72, 2.31, 0;
			       0.97, -2.54, 3.75;
			       1.4, -1.96, 2.5;
			       1.83, -1.39, 1.25;
			       2.26, -0.81, 0;
			       1.94, -5.08, 2.5;
			       2.37, -4.5, 1.25;
			       2.8, -3.92, 0;
			       2.91, -7.62, 1.25;
			       3.34, -7.04, 0;
			       3.88, -10.15, 0];

## Generate 4th order Lagrange basis functions on the reference
## triangle.
lagrange_shape_functions = LagrangeBasisOn3DTria2Args(4);
## Generate a list of area coordinates distributed in the reference
## triangle.
tria_area_coords = AreaCoordsOnTriangle(20);
## Convert the list of area coordinate to global coordinates.
tria_global_coords = ...
AreaToGlobalCoords(tria_area_coords, lagrange_shape_functions,
		   tria_node_global_coord_flat);

figure;
hold on;
Plot3DPointList(tria_global_coords, "r.");
Plot3DPointList(tria_node_global_coord_flat, "bo");
hold off;
