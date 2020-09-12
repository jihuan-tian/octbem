clear all;
ConfigGraphicsToolkit;

addpath("./gmsh_io");
addpath("./stroud");

pkg load fpl;
pkg load msh;

format long;

## Read the mesh.
mesh_filename = "./mesh/sphere-fine-shifted.msh";
problem_domain_mesh = ReadGmshTrias(mesh_filename);
model_sphere_center = [0.2, 0.3, 0.4];
model_sphere_radius = 1.0;

## #############################################################################
## Generate the shape function space for describing the cell geometry, the order
## of which is determined by the mesher.
## #############################################################################
## At present, assume all cells in the mesh have the same geometric
## order, which is only 1st order at the moment.
shape_function_space_order = problem_domain_mesh.cell_geom_orders(1);
## Generate the series of shape functions dependent on the first two components
## of the area coordinate, while the last component is not independent.
shape_function_space = LagrangeBasisOn3DTria2Args(shape_function_space_order);
## Generate the local coordinates for the support nodes of shape functions.
shape_function_support_nodes = AreaCoordsOnTriangle(shape_function_space_order)(:, 1:2);
## Total number of basis functions in the shape function space.
number_of_bases_in_shape_function_space = length(shape_function_space);

## Quadrature norder (number of quadrature points) for the integrals in FEM,
## which will also be used in BEM. N.B. n-point Gauss quadrature is exact for
## polynomials of degree \f$2n - 1\f$.
## Ref: https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature
fem_gauss_quad_norder_2d = ceil((shape_function_space_order + shape_function_space_order + 1) / 2);
## Generate the Gauss-Legendre quadrature rules to perform 2D FEM quadrature.
[fem_gauss_quad_2d_qpts_xi1, fem_gauss_quad_2d_qpts_xi2, fem_gauss_quad_2d_qwts] = triangle_unit_product_set(fem_gauss_quad_norder_2d);
fem_gauss_quad_2d_qpts = [fem_gauss_quad_2d_qpts_xi1', fem_gauss_quad_2d_qpts_xi2'];
fem_gauss_quad_2d_qwts = fem_gauss_quad_2d_qwts';

sphere_area_planar_mesh = 0;
sphere_area_curved_mesh = 0;
for e = 1:problem_domain_mesh.number_of_cells
  cell_node_indices = problem_domain_mesh.mesh_cells(e, :);

  cell_jacobian_functor_planar = @(x) problem_domain_mesh.cell_surface_metrics(e);
  cell_jacobian_functor_curved = @(x) GlobalSurfaceMetricOnS2Tria(x, problem_domain_mesh.mesh_nodes(cell_node_indices, :), shape_function_space, model_sphere_center, model_sphere_radius);
  
  sphere_area_planar_mesh += ApplyQuadRule(@(x) 1, fem_gauss_quad_2d_qpts, fem_gauss_quad_2d_qwts, cell_jacobian_functor_planar) * triangle_unit_volume();
  
  sphere_area_curved_mesh += ApplyQuadRule(@(x) 1, fem_gauss_quad_2d_qpts, fem_gauss_quad_2d_qwts, cell_jacobian_functor_curved) * triangle_unit_volume();
endfor

sphere_area_planar_mesh
sphere_area_curved_mesh
