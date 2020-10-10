## This script tests the quadrature rule proposed in
## Erichsen1996Efficient on high order cell.
clear all;
addpath("./stroud");

## Global coordinates for the nodes of a flat triangle.
mesh_nodes = [0, 0, 5; ...
	      0.43, 0.58, 3.75; ...
	      0.86, 1.15, 2.5; ...
	      1.29, 1.73, 1.25; ...
	      1.72, 2.31, 0; ...
	      0.97, -2.54, 3.75; ...
	      1.4, -1.96, 2.5; ...
	      1.83, -1.39, 1.25; ...
	      2.26, -0.81, 0; ...
	      1.94, -5.08, 2.5; ...
	      2.37, -4.5, 1.25; ...
	      2.8, -3.92, 0; ...
	      2.91, -7.62, 1.25; ...
	      3.34, -7.04, 0; ...
	      3.88, -10.15, 0];
[number_of_nodes, space_dim] = size(mesh_nodes);
mesh_cells = 1:number_of_nodes;
number_of_cells = size(mesh_cells, 1);
mesh_cell_surface_metrics = zeros(number_of_cells, 1);
mesh_cell_normal_vectors = zeros(number_of_cells, space_dim);
mesh_cell_ranges = zeros(number_of_cells, 1);

enable_automatic_quad_order = false;


for m = 1:number_of_cells
  tria_corners = mesh_nodes(mesh_cells(m, GetTriaCornerNodeIndices(TriaGeometryOrder(length(mesh_cells(m, :))))), :);
  mesh_cell_surface_metrics(m) = GlobalSurfaceMetricOn3DFlatTria(tria_corners);
  mesh_cell_normal_vectors(m, :) = SurfaceNormalOn3DFlatTria(tria_corners);
  mesh_cell_ranges(m) = Tria3DRange(tria_corners);
endfor

## Estimation of the mesh size h.
mesh_size_estimate = min(mesh_cell_ranges);

## Sobolev space s1 value.
sobolev_function_space_order = 0;
galerkin_estimate_norm_index = 0;
basis_function_order_for_triangle = 4;
## Generate the shape functions for describing the geometry.
lagrange_shape_functions = LagrangeBasisOn3DTria2Args(basis_function_order_for_triangle);

## Test and basis functions on Kx and Ky (same panel case).
kx_basis_function = lagrange_shape_functions{2};
ky_basis_function = lagrange_shape_functions{7};

## Cell index used for the quadrature.
kx_cell_index = 1;
ky_cell_index = 1;

problem_domain_mesh.mesh_cells = mesh_cells;
problem_domain_mesh.mesh_nodes = mesh_nodes;
problem_domain_mesh.number_of_cells = number_of_cells;
problem_domain_mesh.number_of_nodes = number_of_nodes;

## Calculate the distance between each pair of panels.
panel_distance_matrix = CalcTriaPanelDistanceMatrix(problem_domain_mesh);

## Calculate the neighboring type between each pair of panels.
neighboring_type_matrix = CalcTriaPanelNeighboringTypes(problem_domain_mesh.mesh_cells);

## Construct the functors for normal vector.
nx_functor = @(x) SurfaceNormalOn3DTria(x, mesh_nodes(mesh_cells(kx_cell_index, :), :));
ny_functor = @(x) SurfaceNormalOn3DTria(x, mesh_nodes(mesh_cells(ky_cell_index, :), :));

## Verify the inner product between the nx_functor and the constant value nx.
for xi = 0:0.05:1
 for eta = 0:0.05:1
   inner_prod(nx_functor([xi, eta]), mesh_cell_normal_vectors(1, :))
 endfor
endfor

## Construct the functors for surface metric.
Jx_functor = @(x) GlobalSurfaceMetricOn3DTria(x, mesh_nodes(mesh_cells(kx_cell_index, :), :));
Jy_functor = @(x) GlobalSurfaceMetricOn3DTria(x, mesh_nodes(mesh_cells(ky_cell_index, :), :));

## SLP
kernel_singularity_order = 1;
laplace_slp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceSLPKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(1, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(1), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

laplace_slp_kernel_3d_integral = ErichsenQuadRule(@LaplaceSLPKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, nx_functor, ny_functor, Jx_functor, Jy_functor, mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## DLP
kernel_singularity_order = 2;
laplace_dlp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceDLPKernel3DFlat, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(1, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(1), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

laplace_dlp_kernel_3d_integral = ErichsenQuadRule(@LaplaceDLPKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, nx_functor, ny_functor, Jx_functor, Jy_functor, mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## ADLP
kernel_singularity_order = 2;
laplace_adlp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceDLPAdjointKernel3DFlat, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(1, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(1), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

laplace_adlp_kernel_3d_integral = ErichsenQuadRule(@LaplaceDLPAdjointKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, nx_functor, ny_functor, Jx_functor, Jy_functor, mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## Hyper-sing
kernel_singularity_order = 3;
laplace_hyper_sing_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceHyperSingKernel3DFlat, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(1, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(1), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

laplace_hyper_sing_kernel_3d_integral = ErichsenQuadRule(@LaplaceHyperSingKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, nx_functor, ny_functor, Jx_functor, Jy_functor, mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)
