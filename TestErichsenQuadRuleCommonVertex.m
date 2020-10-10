clear all;
addpath("./stroud");

mesh_nodes = [-8, 2, 0;
	      -5.085184886164775, -4.356393710357985, -2.557625847099054;
	      -2.282644513289576, 1.271453055760844, 2;
	      3.332081929919845, 8.10494765534936, 1.475185518874472;
	      4.911995445557295, 0.721841246399158, 2.659304210615577];
mesh_cells = [3, 1, 2;
	      4, 3, 5];
[number_of_nodes, space_dim] = size(mesh_nodes);
number_of_cells = size(mesh_cells, 1);
mesh_cell_surface_metrics = zeros(number_of_cells, 1);
mesh_cell_normal_vectors = zeros(number_of_cells, space_dim);
mesh_cell_ranges = zeros(number_of_cells, 1);

for m = 1:number_of_cells
  mesh_cell_surface_metrics(m) = GlobalSurfaceMetricOn3DFlatTria(mesh_nodes(mesh_cells(m, :), :));
  mesh_cell_normal_vectors(m, :) = SurfaceNormalOn3DFlatTria(mesh_nodes(mesh_cells(m, :), :));
  mesh_cell_ranges(m) = Tria3DRange(mesh_nodes(mesh_cells(m, :), :));
endfor

## Estimation of the mesh size h.
mesh_size_estimate = min(mesh_cell_ranges);

## Sobolev space s1 value.
sobolev_function_space_order = 0;
galerkin_estimate_norm_index = 0;
basis_function_order_for_triangle = 1;
## Generate shape functions for describing the geometry.
lagrange_shape_functions = LagrangeBasisOn3DTria2Args(basis_function_order_for_triangle);

## Test and basis functions on Kx and Ky.
kx_basis_function = lagrange_shape_functions{1};
ky_basis_function = lagrange_shape_functions{1};

## Cell index used for the quadrature.
kx_cell_index = 1;
ky_cell_index = 2;

problem_domain_mesh.mesh_cells = mesh_cells;
problem_domain_mesh.mesh_nodes = mesh_nodes;
problem_domain_mesh.number_of_cells = number_of_cells;
problem_domain_mesh.number_of_nodes = number_of_nodes;

## Calculate the distance between each pair of panels.
panel_distance_matrix = CalcTriaPanelDistanceMatrix(problem_domain_mesh);

## Calculate the neighboring type between each pair of panels.
neighboring_type_matrix = CalcTriaPanelNeighboringTypes(problem_domain_mesh.mesh_cells);

## SLP
kernel_singularity_order = 1;
laplace_slp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceSLPKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## DLP
kernel_singularity_order = 2;
laplace_dlp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceDLPKernel3DFlat, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## ADLP
kernel_singularity_order = 2;
laplace_adlp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceDLPAdjointKernel3DFlat, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## Hyper-sing
kernel_singularity_order = 3;
laplace_hyper_sing_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceHyperSingKernel3DFlat, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, panel_distance_matrix, neighboring_type_matrix, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)
