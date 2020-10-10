ConfigGraphicsToolkit;

## ##############################################################################
## Evaluate each Lagrange basis functions with 3 arguments on its own supporting
## node.
## ##############################################################################
max_basis_function_order = 3;
for basis_function_order = 0:max_basis_function_order
  lagrange_basis_functions = LagrangeBasisOn3DTria3Args(basis_function_order);
  tria_area_coords = AreaCoordsOnTriangle(basis_function_order);
  number_of_lagrange_basis_functions = NumberOfNodesOnTriangle(basis_function_order);
  lagrange_basis_function_values = zeros(number_of_lagrange_basis_functions, number_of_lagrange_basis_functions);

  for m = 1:number_of_lagrange_basis_functions
    lagrange_basis_function_values(:, m) = lagrange_basis_functions{m}(tria_area_coords);
  endfor

  NewFigure;
  colormap("copper");
  imagesc(lagrange_basis_function_values);
  colorbar;
  view(2);
  axis equal;
  if (basis_function_order > 0)
    xlim([1, number_of_lagrange_basis_functions]);
    ylim([1, number_of_lagrange_basis_functions]);
  endif
  xlabel("Lagrange basis function No.");
  ylabel("Node No. in reference triangle");
  title(cstrcat("Evaluation of Lagrange basis functions with order ", num2str(basis_function_order)));
endfor

## ##############################################################################
## Evaluate each Lagrange basis functions with 2 arguments on its own supporting
## node.
## ##############################################################################
max_basis_function_order = 3;
for basis_function_order = 1:max_basis_function_order
  lagrange_basis_functions = LagrangeBasisOn3DTria2Args(basis_function_order);
  tria_area_coords = AreaCoordsOnTriangle(basis_function_order);
  number_of_lagrange_basis_functions = NumberOfNodesOnTriangle(basis_function_order);
  lagrange_basis_function_values = zeros(number_of_lagrange_basis_functions, number_of_lagrange_basis_functions);

  for m = 1:number_of_lagrange_basis_functions
    lagrange_basis_function_values(:, m) = lagrange_basis_functions{m}(tria_area_coords(:,1:2));
  endfor

  NewFigure;
  colormap("copper");
  imagesc(lagrange_basis_function_values);
  colorbar;
  view(2);
  axis equal;
  if (basis_function_order > 0)
    xlim([1, number_of_lagrange_basis_functions]);
    ylim([1, number_of_lagrange_basis_functions]);
  endif
  xlabel("Lagrange basis function No.");
  ylabel("Node No. in reference triangle");
  title(cstrcat("Evaluation of Lagrange basis functions with order ", num2str(basis_function_order)));
endfor

## #################################################
## Calculate the Jacobi matrix and its determinant.
## #################################################
basis_function_order_for_triangle = 4;

## Global coordinates for the nodes of a flat triangle.
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

## Global coordinates for the nodes of a curved triangle.
tria_node_global_coord_curved = [-4.01406, 9.94323, 4.5;
				 -3.65318, 6.16233, 3.8;
				 -2.29768, 2.279, 3;
				 -0.282742, -0.981539, 3.8;
				 2.8038, -4.13172, 4.5;
				 -6.69391, 6.71186, 3.8;
				 -5.41167, 2.64535, 3;
				 -3.14029, -1.5677, 3;
				 -0.575822, -4.82824, 3.8;
				 -8.41576, 3.48796, 3;
				 -6.40083, -0.725094, 2;
				 -3.87299, -4.82824, 3;
				 -9.58809, 0.00761196, 2;
				 -7.13353, -4.24207, 1;
				 -10.0638, -3.29188, 1];

## Evalute the Jacobi matrix at each node in the reference cell.
tria_area_coords = AreaCoordsOnTriangle(basis_function_order_for_triangle);
shape_function_jacobi_matrix_3args = ShapeFunctionJacobiOn3DTria3Args(basis_function_order_for_triangle);
shape_function_jacobi_matrix_2args = ShapeFunctionJacobiOn3DTria2Args(basis_function_order_for_triangle);

## Jacobi matrix for the flat triangle
jacobi_matrix_3args_for_flat_tria_yz = GlobalJacobiMatrix(tria_area_coords, [tria_node_global_coord_flat(:,2), tria_node_global_coord_flat(:,3)], shape_function_jacobi_matrix_3args);
jacobi_matrix_3args_for_flat_tria_zx = GlobalJacobiMatrix(tria_area_coords, [tria_node_global_coord_flat(:,3), tria_node_global_coord_flat(:,1)], shape_function_jacobi_matrix_3args);
jacobi_matrix_3args_for_flat_tria_xy = GlobalJacobiMatrix(tria_area_coords, [tria_node_global_coord_flat(:,1), tria_node_global_coord_flat(:,2)], shape_function_jacobi_matrix_3args);

jacobi_matrix_2args_for_flat_tria_yz = GlobalJacobiMatrix(tria_area_coords(:,1:2), [tria_node_global_coord_flat(:,2), tria_node_global_coord_flat(:,3)], shape_function_jacobi_matrix_2args);
jacobi_matrix_2args_for_flat_tria_zx = GlobalJacobiMatrix(tria_area_coords(:,1:2), [tria_node_global_coord_flat(:,3), tria_node_global_coord_flat(:,1)], shape_function_jacobi_matrix_2args);
jacobi_matrix_2args_for_flat_tria_xy = GlobalJacobiMatrix(tria_area_coords(:,1:2), [tria_node_global_coord_flat(:,1), tria_node_global_coord_flat(:,2)], shape_function_jacobi_matrix_2args);

## Jacobi matrix for the curved triangle
jacobi_matrix_3args_for_curved_tria_yz = GlobalJacobiMatrix(tria_area_coords, [tria_node_global_coord_curved(:,2), tria_node_global_coord_curved(:,3)], shape_function_jacobi_matrix_3args);
jacobi_matrix_3args_for_curved_tria_zx = GlobalJacobiMatrix(tria_area_coords, [tria_node_global_coord_curved(:,3), tria_node_global_coord_curved(:,1)], shape_function_jacobi_matrix_3args);
jacobi_matrix_3args_for_curved_tria_xy = GlobalJacobiMatrix(tria_area_coords, [tria_node_global_coord_curved(:,1), tria_node_global_coord_curved(:,2)], shape_function_jacobi_matrix_3args);

jacobi_matrix_2args_for_curved_tria_yz = GlobalJacobiMatrix(tria_area_coords(:,1:2), [tria_node_global_coord_curved(:,2), tria_node_global_coord_curved(:,3)], shape_function_jacobi_matrix_2args);
jacobi_matrix_2args_for_curved_tria_zx = GlobalJacobiMatrix(tria_area_coords(:,1:2), [tria_node_global_coord_curved(:,3), tria_node_global_coord_curved(:,1)], shape_function_jacobi_matrix_2args);
jacobi_matrix_2args_for_curved_tria_xy = GlobalJacobiMatrix(tria_area_coords(:,1:2), [tria_node_global_coord_curved(:,1), tria_node_global_coord_curved(:,2)], shape_function_jacobi_matrix_2args);

## Jacobi det for the flat triangle
jacobi_det_3args_for_flat_tria_yz = GlobalJacobiDet(tria_area_coords, [tria_node_global_coord_flat(:,2), tria_node_global_coord_flat(:,3)], shape_function_jacobi_matrix_3args);
jacobi_det_3args_for_flat_tria_zx = GlobalJacobiDet(tria_area_coords, [tria_node_global_coord_flat(:,3), tria_node_global_coord_flat(:,1)], shape_function_jacobi_matrix_3args);
jacobi_det_3args_for_flat_tria_xy = GlobalJacobiDet(tria_area_coords, [tria_node_global_coord_flat(:,1), tria_node_global_coord_flat(:,2)], shape_function_jacobi_matrix_3args);

jacobi_det_2args_for_flat_tria_yz = GlobalJacobiDet(tria_area_coords(:,1:2), [tria_node_global_coord_flat(:,2), tria_node_global_coord_flat(:,3)], shape_function_jacobi_matrix_2args);
jacobi_det_2args_for_flat_tria_zx = GlobalJacobiDet(tria_area_coords(:,1:2), [tria_node_global_coord_flat(:,3), tria_node_global_coord_flat(:,1)], shape_function_jacobi_matrix_2args);
jacobi_det_2args_for_flat_tria_xy = GlobalJacobiDet(tria_area_coords(:,1:2), [tria_node_global_coord_flat(:,1), tria_node_global_coord_flat(:,2)], shape_function_jacobi_matrix_2args);

## Jacobi det for the curved triangle
jacobi_det_3args_for_curved_tria_yz = GlobalJacobiDet(tria_area_coords, [tria_node_global_coord_curved(:,2), tria_node_global_coord_curved(:,3)], shape_function_jacobi_matrix_3args);
jacobi_det_3args_for_curved_tria_zx = GlobalJacobiDet(tria_area_coords, [tria_node_global_coord_curved(:,3), tria_node_global_coord_curved(:,1)], shape_function_jacobi_matrix_3args);
jacobi_det_3args_for_curved_tria_xy = GlobalJacobiDet(tria_area_coords, [tria_node_global_coord_curved(:,1), tria_node_global_coord_curved(:,2)], shape_function_jacobi_matrix_3args);

jacobi_det_2args_for_curved_tria_yz = GlobalJacobiDet(tria_area_coords(:,1:2), [tria_node_global_coord_curved(:,2), tria_node_global_coord_curved(:,3)], shape_function_jacobi_matrix_2args);
jacobi_det_2args_for_curved_tria_zx = GlobalJacobiDet(tria_area_coords(:,1:2), [tria_node_global_coord_curved(:,3), tria_node_global_coord_curved(:,1)], shape_function_jacobi_matrix_2args);
jacobi_det_2args_for_curved_tria_xy = GlobalJacobiDet(tria_area_coords(:,1:2), [tria_node_global_coord_curved(:,1), tria_node_global_coord_curved(:,2)], shape_function_jacobi_matrix_2args);

## Calculate the surface metric for the flat triangle: its value should be 2
## times of the triangle area: 34.28272508151879.
surface_metric_3args_for_flat_tria = GlobalSurfaceMetricOn3DTria(tria_area_coords, tria_node_global_coord_flat);
surface_metric_2args_for_flat_tria = GlobalSurfaceMetricOn3DTria(tria_area_coords(:,1:2), tria_node_global_coord_flat);

## Calculate the surface metric for the curved triangle.
surface_metric_3args_for_curved_tria = GlobalSurfaceMetricOn3DTria(tria_area_coords, tria_node_global_coord_curved);
surface_metric_2args_for_curved_tria = GlobalSurfaceMetricOn3DTria(tria_area_coords(:,1:2), tria_node_global_coord_curved);

## Calculate the surface normals for the flat triangle.
surface_normals_3args_for_flat_tria = SurfaceNormalOn3DTria(tria_area_coords, tria_node_global_coord_flat);
surface_normals_2args_for_flat_tria = SurfaceNormalOn3DTria(tria_area_coords(:,1:2), tria_node_global_coord_flat);

## Calculate the surface normals for the curved triangle.
surface_normals_3args_for_curved_tria = SurfaceNormalOn3DTria(tria_area_coords, tria_node_global_coord_curved);
surface_normals_2args_for_curved_tria = SurfaceNormalOn3DTria(tria_area_coords(:,1:2), tria_node_global_coord_curved);

## Visualize the area coordinates and their corresponding normal vectors.
NewFigure;
hold on;
Plot3DPointList(tria_node_global_coord_flat, "r-");
Plot3DAffineVectors(tria_node_global_coord_flat, surface_normals_3args_for_flat_tria);

NewFigure;
hold on;
Plot3DAffineVectors(tria_node_global_coord_curved, surface_normals_3args_for_curved_tria);
Plot3DPointList(tria_node_global_coord_curved, "bo");
