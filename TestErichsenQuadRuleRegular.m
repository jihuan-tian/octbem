clear all;
addpath("./stroud");

mesh_nodes = [-8, 2, 0;
	      -5.085184886164775, -4.356393710357985, -2.557625847099054;
	      -2.282644513289576, 1.271453055760844, 2;
	      -27.49303083484554, 17.091128129664433, 13.508059530450183;
	      -34.64806071646211, 12.124069372654965, 11.657881608461125;
	      -26.01454406448393, 5.889056159725726, 17.688754663650442];
mesh_cells = [3, 1, 2;
	      6, 4, 5];
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

## SLP
kernel_singularity_order = 1;
laplace_slp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceSLPKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## DLP
kernel_singularity_order = 2;
laplace_dlp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceDLPKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## ADLP
kernel_singularity_order = 2;
laplace_adlp_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceDLPAdjointKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)

## Hyper-sing
kernel_singularity_order = 3;
laplace_hyper_sing_kernel_3d_integral = ErichsenQuadRuleFlat(@LaplaceHyperSingKernel3D, kernel_singularity_order, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, kx_cell_index, ky_cell_index, mesh_cells, mesh_nodes, mesh_cell_normal_vectors(1, :), mesh_cell_normal_vectors(2, :), mesh_cell_surface_metrics(1), mesh_cell_surface_metrics(2), mesh_size_estimate, sobolev_function_space_order, basis_function_order_for_triangle, galerkin_estimate_norm_index)
