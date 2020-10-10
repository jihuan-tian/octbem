ConfigGraphicsToolkit;

## ##############################################################################
## Evaluate each Lagrange basis functions on its own supporting node.
## ##############################################################################
max_basis_function_order = 3;
for basis_function_order = 1:max_basis_function_order
  lagrange_basis_functions = LagrangeBasisOn3DQuad(basis_function_order);
  quad_area_coords = AreaCoordsOnQuad(basis_function_order);
  number_of_lagrange_basis_functions = NumberOfNodesOnQuadrangle(basis_function_order);
  lagrange_basis_function_values = zeros(number_of_lagrange_basis_functions, number_of_lagrange_basis_functions);

  for m = 1:number_of_lagrange_basis_functions
    lagrange_basis_function_values(:, m) = lagrange_basis_functions{m}(quad_area_coords);
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
  ylabel("Node No. in reference quadrangle");
  title(cstrcat("Evaluation of Lagrange basis functions with order ", num2str(basis_function_order)));
endfor

## ## Plot the Lagrange function distribution.
## lagrange_basis_function_values = zeros(100, 100);
## xi = linspace(0, 1, 100);
## eta = xi;
## for k = 1:length(lagrange_basis_functions)
##   for m = 1:100
##     for n = 1:100
##       lagrange_basis_function_values(m, n) = lagrange_basis_functions{k}([xi(m), eta(n)]);
##     endfor
##   endfor
##   NewFigure;
##   colormap("copper");
##   imagesc(lagrange_basis_function_values);
##   colorbar;
##   view(2);
##   axis equal;
## endfor

## #################################################
## Calculate the Jacobi matrix and its determinant.
## #################################################
basis_function_order_for_quadrangle = 1;

## Global coordinates for the nodes of a flat quadrangle.
quad_node_global_coord_flat = [-8.56704360605347, 5.10054617887335, 2.94480170938683; ...
			       -5.48877539622046, 6.63967635136584, 3.83341892879241; ...
			       4.50996329004166, -7.73179943406246, -4.46395648457663; ...
			       -4.51627282749619, -5.11395296568278, -2.95254212136043];

## N.B. Swap the 3rd and 4th nodes in the cell to be consistent with
## the zigzag traversing convention.
temp_node = quad_node_global_coord_flat(4, :);
quad_node_global_coord_flat(4, :) = quad_node_global_coord_flat(3, :);
quad_node_global_coord_flat(3, :) = temp_node;

## Evalute the Jacobi matrix at each node in the reference cell.
quad_area_coords = AreaCoordsOnQuad(basis_function_order_for_quadrangle);
shape_function_jacobi_matrix = ShapeFunctionJacobiOn3DQuad(basis_function_order_for_quadrangle);

## Jacobi matrix for the flat quadrangle.
jacobi_matrix_for_flat_quad_yz = GlobalJacobiMatrix(quad_area_coords, [quad_node_global_coord_flat(:,2), quad_node_global_coord_flat(:,3)], shape_function_jacobi_matrix);
jacobi_matrix_for_flat_quad_zx = GlobalJacobiMatrix(quad_area_coords, [quad_node_global_coord_flat(:,3), quad_node_global_coord_flat(:,1)], shape_function_jacobi_matrix);
jacobi_matrix_for_flat_quad_xy = GlobalJacobiMatrix(quad_area_coords, [quad_node_global_coord_flat(:,1), quad_node_global_coord_flat(:,2)], shape_function_jacobi_matrix);

## Jacobi det for the flat quadrangle.
jacobi_det_for_flat_quad_yz = GlobalJacobiDet(quad_area_coords, [quad_node_global_coord_flat(:,2), quad_node_global_coord_flat(:,3)], shape_function_jacobi_matrix);
jacobi_det_for_flat_quad_zx = GlobalJacobiDet(quad_area_coords, [quad_node_global_coord_flat(:,3), quad_node_global_coord_flat(:,1)], shape_function_jacobi_matrix);
jacobi_det_for_flat_quad_xy = GlobalJacobiDet(quad_area_coords, [quad_node_global_coord_flat(:,1), quad_node_global_coord_flat(:,2)], shape_function_jacobi_matrix);

## Calculate the surface metric for the flat quadrangle.
surface_metric_for_flat_quad = GlobalSurfaceMetricOn3DQuad(quad_area_coords, quad_node_global_coord_flat);

## Calculate the surface normals for the flat quadrangle.
## The normal vector from GeoGebra is: [0.000000000000000e+00, 5.000000000005073e-01, -8.660254037841457e-01], which is just the negation of the calculation results here.
surface_normals_for_flat_quad = SurfaceNormalOn3DQuad(quad_area_coords, quad_node_global_coord_flat);

## Calculate the area of the flat quadrangle using numerical quadrature.
norder_for_quadrature = 3;
[qpts_1d, qwts_1d] = GaussLegendreRule(norder_for_quadrature);
## Adjust the Gauss-Legendre rule from [-1, 1] to [0, 1].
[qpts_1d, qwts_1d] = rule_adjust(-1, 1, 0, 1, norder_for_quadrature, qpts_1d, qwts_1d);
qpts_2d = tensor_prod2(qpts_1d, qpts_1d);
qwts_2d = tensor_prod2(qwts_1d, qwts_1d);
quad_area = ApplyQuadRule(@(x) GlobalSurfaceMetricOn3DQuad(x, quad_node_global_coord_flat), qpts_2d, qwts_2d)

## Visualize the area coordinates and their corresponding normal vectors.
NewFigure;
hold on;
Plot3DPointList(quad_node_global_coord_flat, "r-");
Plot3DAffineVectors(quad_node_global_coord_flat, surface_normals_for_flat_quad);
xlabel("X");
ylabel("Y");
zlabel("Z");
axis equal;
