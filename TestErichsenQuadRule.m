## This script tests the quadrature rule proposed in Erichsen1996Efficient.

clear all;
addpath("./stroud");

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

enable_automatic_quad_order = false;
tria_corners = tria_node_global_coord_flat([1, 5, 15], :);
tria_circumcenter = Tria3DCircumcenter(tria_corners);
tria_range = Tria3DRange(tria_corners);
tria_surface_metric = GlobalSurfaceMetricOn3DFlatTria(tria_corners);
tria_normal_vector = SurfaceNormalOn3DFlatTria(tria_corners);
## Sobolev space s1 value.
function_space_order = 0;
if (enable_automatic_quad_order)
  norder_for_eta = ceil((basis_function_order_for_triangle - function_space_order) / 2 * abs(log(tria_range)));
  norder_for_omega = ceil(basis_function_order_for_triangle - function_space_order);
else
  ## The quadrature order is given on the command line.
  norder_for_eta = 6;
  norder_for_omega = 6;
endif

## Generate the shape functions.
lagrange_shape_functions = LagrangeBasisOn3DTria2Args(basis_function_order_for_triangle);

## Test and basis functions on Kx and Ky (same panel case).
kx_basis_function = lagrange_shape_functions{2};
ky_basis_function = lagrange_shape_functions{7};

## ############################################################################
## Calculate the integration for all kinds of kernels with different quadrature
## orders.
## ############################################################################
laplace_slp_kernel_3d_integral = ErichsenQuadSamePanel(@LaplaceSLPKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, @SurfaceNormalOn3DFlatTria, @SurfaceNormalOn3DFlatTria);
laplace_dlp_kernel_3d_integral = ErichsenQuadSamePanel(@LaplaceDLPKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, @SurfaceNormalOn3DTria, @SurfaceNormalOn3DTria);
laplace_adlp_kernel_3d_integral = ErichsenQuadSamePanel(@LaplaceDLPAdjointKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, @SurfaceNormalOn3DTria, @SurfaceNormalOn3DTria);
laplace_hyper_kernel_3d_integral = ErichsenQuadSamePanel(@LaplaceHyperSingKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, @SurfaceNormalOn3DTria, @SurfaceNormalOn3DTria);

laplace_slp_kernel_3d_integral = ErichsenQuadSamePanelFlat(@LaplaceSLPKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, tria_normal_vector, tria_normal_vector, tria_surface_metric, tria_surface_metric);
laplace_dlp_kernel_3d_integral = ErichsenQuadSamePanelFlat(@LaplaceDLPKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, tria_normal_vector, tria_normal_vector, tria_surface_metric, tria_surface_metric);
laplace_adlp_kernel_3d_integral = ErichsenQuadSamePanelFlat(@LaplaceDLPAdjointKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, tria_normal_vector, tria_normal_vector, tria_surface_metric, tria_surface_metric);
laplace_hyper_kernel_3d_integral = ErichsenQuadSamePanelFlat(@LaplaceHyperSingKernel3D, norder_for_eta, norder_for_omega, kx_basis_function, ky_basis_function, lagrange_shape_functions, lagrange_shape_functions, tria_node_global_coord_flat, tria_node_global_coord_flat, tria_normal_vector, tria_normal_vector, tria_surface_metric, tria_surface_metric);
