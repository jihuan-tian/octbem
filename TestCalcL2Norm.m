clear all;
ConfigGraphicsToolkit;

addpath("./gmsh_io");
addpath("./stroud");

pkg load fpl;
pkg load msh;

mesh_filename = "./mesh/sphere-coarse.msh";
problem_domain_mesh = ReadGmshTrias(mesh_filename);
ansatz_function_space_order = problem_domain_mesh.cell_geom_orders(1);
ansatz_function_space = LagrangeBasisOn3DTria2Args(ansatz_function_space_order);

fem_gauss_quad_norder_2d = ceil((ansatz_function_space_order + ansatz_function_space_order + 1) / 2);
## Generate the Gauss-Legendre quadrature rules to perform 2D FEM quadrature.
[fem_gauss_quad_2d_qpts_xi1, fem_gauss_quad_2d_qpts_xi2, fem_gauss_quad_2d_qwts] = triangle_unit_product_set(fem_gauss_quad_norder_2d);
fem_gauss_quad_2d_qpts = [fem_gauss_quad_2d_qpts_xi1', fem_gauss_quad_2d_qpts_xi2'];
fem_gauss_quad_2d_qwts = fem_gauss_quad_2d_qwts';

CalcL2Norm(ones(problem_domain_mesh.number_of_nodes), ansatz_function_space, 1:problem_domain_mesh.number_of_cells, problem_domain_mesh, fem_gauss_quad_norder_2d)

mesh_filename = "./mesh/sphere.msh";
problem_domain_mesh = ReadGmshTrias(mesh_filename);
ansatz_function_space_order = problem_domain_mesh.cell_geom_orders(1);
ansatz_function_space = LagrangeBasisOn3DTria2Args(ansatz_function_space_order);

fem_gauss_quad_norder_2d = ceil((ansatz_function_space_order + ansatz_function_space_order + 1) / 2);
## Generate the Gauss-Legendre quadrature rules to perform 2D FEM quadrature.
[fem_gauss_quad_2d_qpts_xi1, fem_gauss_quad_2d_qpts_xi2, fem_gauss_quad_2d_qwts] = triangle_unit_product_set(fem_gauss_quad_norder_2d);
fem_gauss_quad_2d_qpts = [fem_gauss_quad_2d_qpts_xi1', fem_gauss_quad_2d_qpts_xi2'];
fem_gauss_quad_2d_qwts = fem_gauss_quad_2d_qwts';

CalcL2Norm(ones(problem_domain_mesh.number_of_nodes), ansatz_function_space, 1:problem_domain_mesh.number_of_cells, problem_domain_mesh, fem_gauss_quad_norder_2d)

mesh_filename = "./mesh/sphere-fine.msh";
problem_domain_mesh = ReadGmshTrias(mesh_filename);
ansatz_function_space_order = problem_domain_mesh.cell_geom_orders(1);
ansatz_function_space = LagrangeBasisOn3DTria2Args(ansatz_function_space_order);

fem_gauss_quad_norder_2d = ceil((ansatz_function_space_order + ansatz_function_space_order + 1) / 2);
## Generate the Gauss-Legendre quadrature rules to perform 2D FEM quadrature.
[fem_gauss_quad_2d_qpts_xi1, fem_gauss_quad_2d_qpts_xi2, fem_gauss_quad_2d_qwts] = triangle_unit_product_set(fem_gauss_quad_norder_2d);
fem_gauss_quad_2d_qpts = [fem_gauss_quad_2d_qpts_xi1', fem_gauss_quad_2d_qpts_xi2'];
fem_gauss_quad_2d_qwts = fem_gauss_quad_2d_qwts';

CalcL2Norm(ones(problem_domain_mesh.number_of_nodes), ansatz_function_space, 1:problem_domain_mesh.number_of_cells, problem_domain_mesh, fem_gauss_quad_norder_2d)
