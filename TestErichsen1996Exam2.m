clear all;
ConfigGraphicsToolkit;

addpath("./gmsh_io");
addpath("./stroud");

pkg load fpl;
pkg load msh;

## Read the mesh.
mesh_filename = "./mesh/sphere-coarse.msh";
problem_domain_mesh = ReadGmshTrias(mesh_filename);
## Plot the mesh.
figure;
triplot3d(problem_domain_mesh.mesh_cells, problem_domain_mesh.mesh_nodes, "color", "r", "linestyle", "-", "marker", "o");
hold on;
Plot3DAffineVectors(problem_domain_mesh.cell_cogs, problem_domain_mesh.cell_normal_vectors / 4);
grid off;
axis equal;
axis off;
hold off;

## Calculate the distance of the mesh nodes to the origin.
fprintf(stdout(), "Distance of mesh nodes to the origin:\n");
for m = 1:problem_domain_mesh.number_of_nodes
  nodes_to_origin_dists(m) = norm(problem_domain_mesh.mesh_nodes(m, :));
endfor

## Calculate the distance of the cell's center of gravity to the origin.
fprintf(stdout(), "Distance of cell's center of gravity to the origin:\n");
for m = 1:problem_domain_mesh.number_of_cells
  cogs_to_origin_dists(m) = norm(problem_domain_mesh.cell_cogs(m, :));
endfor

## Calculate the distance of the end of surface normal vectors to the
## origin, in order to verify that they actuall point outwards, i.e. the calculated distance should be larger than the cell's center of gravity to the origin.
fprintf(stdout(), "Distance of cell surface normal vector end points to the origin:\n");
cell_normal_vector_endpoints = (problem_domain_mesh.cell_cogs + problem_domain_mesh.cell_normal_vectors);
for m = 1:problem_domain_mesh.number_of_cells
  cell_normal_to_origin_dists(m) = norm(cell_normal_vector_endpoints(m, :));
endfor

if (!isempty(find((cell_normal_to_origin_dists > cogs_to_origin_dists) == 0)))
  error("There is a inward cell's surface normal vector!");
endif

## Location of the unit Dirac point source.
x0 = [0.25, 0.25, 0.25];
## The exact solution:
u_star = zeros(problem_domain_mesh.number_of_nodes, 1);
## Boundary condition on the domain surface (assigned to each node):
## \f$\frac{\pdiff u}{\pdiff n}\f$, which is the normal derivative.
neumann_bc = zeros(problem_domain_mesh.number_of_nodes, 1);
## Assign boundary condition and generate the analytical solution on the nodes.
for m = 1:problem_domain_mesh.number_of_nodes
  ## \f$x - x_0\f$
  diff_vector = problem_domain_mesh.mesh_nodes(m, :) - x0;
  
  ## Analytical solution of the problem.
  u_star(m) = 1.0 / 4.0 / pi / norm(diff_vector);

  ## N.B. The unit normal vector on the spherical surface should point into the
  ## sphere. Because \f$x\f$ lies on the unit spherical surface, \f$-x\f$ is
  ## just the said unit normal vector located at x.
  neumann_bc(m) = inner_prod(problem_domain_mesh.mesh_nodes(m, :), diff_vector) / 4.0 / pi / norm(diff_vector)^3;
endfor

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

## #############################################################################
## Generate the test and ansatz(trial) function spaces: the Lagrange functions
## are used.
## #############################################################################
## At present, the code for generating high order nodes from the lower nodes
## (directly obtained from the mesher) is not implemented. Therefore, the test
## function space order is the same as the shape function space order, which is
## constructed directly from the mesher's output.
test_function_space_order = shape_function_space_order;
test_function_space = LagrangeBasisOn3DTria2Args(test_function_space_order);
## Generate the support nodes associated with basis functions in the test
## function space, keeping the same sequence as the basis functions.
test_function_support_nodes = AreaCoordsOnTriangle(test_function_space_order)(:, 1:2);
## Total number of basis functions in the test function space.
number_of_bases_in_test_function_space = length(test_function_space);

## In Galerkin method, the ansatz(trial) function space is selected to be the
## same as test function space.
ansatz_function_space_order = test_function_space_order;
ansatz_function_space = test_function_space;
ansatz_function_support_nodes = test_function_support_nodes;
number_of_bases_in_ansatz_function_space = length(ansatz_function_space);

## Quadrature norder (number of quadrature points) for the integrals in FEM,
## which will also be used in BEM. N.B. n-point Gauss quadrature is exact for
## polynomials of degree \f$2n - 1\f$.
## Ref: https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature
fem_gauss_quad_norder_2d = ceil((test_function_space_order + ansatz_function_space_order + 1) / 2);
## Generate the Gauss-Legendre quadrature rules to perform 2D FEM quadrature.
[fem_gauss_quad_2d_qpts_xi1, fem_gauss_quad_2d_qpts_xi2, fem_gauss_quad_2d_qwts] = triangle_unit_product_set(fem_gauss_quad_norder_2d);
fem_gauss_quad_2d_qpts = [fem_gauss_quad_2d_qpts_xi1', fem_gauss_quad_2d_qpts_xi2'];
fem_gauss_quad_2d_qwts = fem_gauss_quad_2d_qwts';

stiffness_matrix = zeros(problem_domain_mesh.number_of_nodes, problem_domain_mesh.number_of_nodes);
rhs_matrix = zeros(problem_domain_mesh.number_of_nodes, problem_domain_mesh.number_of_nodes);

sobolev_function_space_order = 0;
galerkin_estimate_norm_index = 0;

## Calculate \f$(v, \frac{1}{2}u)\f$.
for e = 1:problem_domain_mesh.number_of_cells
  local_matrix = zeros(number_of_bases_in_test_function_space, number_of_bases_in_ansatz_function_space);
  for i = 1:number_of_bases_in_test_function_space
    for j = 1:number_of_bases_in_ansatz_function_space
      ## Get the cell Jacobian.
      cell_jacobian_functor = @(x) problem_domain_mesh.cell_surface_metrics(e);
      local_matrix(i, j) = ApplyQuadRule(@(area_coord) 0.5 * test_function_space{i}(area_coord) * ansatz_function_space{j}(area_coord), fem_gauss_quad_2d_qpts, fem_gauss_quad_2d_qwts, cell_jacobian_functor) * triangle_unit_volume();
      
      ## Assemble the local matrix to the global matrix.
      stiffness_matrix(problem_domain_mesh.mesh_cells(e, i), problem_domain_mesh.mesh_cells(e, j)) += local_matrix(i, j);
    endfor
  endfor
endfor

## Calculate \f$v, Ku\f$ and \f$V(\psi)\f$.
## N.B. The default normal vector of the cell from the mesher points outward.
## Its direction should be reversed.
progress_counter = 0;
progress_total_steps = problem_domain_mesh.number_of_cells * problem_domain_mesh.number_of_cells;

## Calculate the distance between each pair of panels.
panel_distance_matrix = CalcPanelDistanceMatrix(problem_domain_mesh);

## Calculate the neighboring type between each pair of panels.
neighboring_type_matrix = CalcPanelNeighboringTypes(problem_domain_mesh.mesh_cells);

## Iterate over each field cell.
for e = 1:problem_domain_mesh.number_of_cells
  ## Iterate over each source cell.
  for f = 1:problem_domain_mesh.number_of_cells
    ## Iterate over each test function, which is associated with the
    ## field cell.
    for i = 1:number_of_bases_in_test_function_space
      ## Iterate over each ansatz function, which is associated with
      ## the source cell.
      for j = 1:number_of_bases_in_ansatz_function_space
	stiffness_matrix(problem_domain_mesh.mesh_cells(e, i), problem_domain_mesh.mesh_cells(f, j)) +=	ErichsenQuadRuleFlat(@LaplaceDLPKernel3DFlat, 2, test_function_space{i}, ansatz_function_space{j}, shape_function_space, shape_function_space, e, f, panel_distance_matrix, neighboring_type_matrix, problem_domain_mesh.mesh_cells, problem_domain_mesh.mesh_nodes, -problem_domain_mesh.cell_normal_vectors(e, :), -problem_domain_mesh.cell_normal_vectors(f, :), problem_domain_mesh.cell_surface_metrics(e), problem_domain_mesh.cell_surface_metrics(f), problem_domain_mesh.min_cell_range, sobolev_function_space_order, test_function_space_order, galerkin_estimate_norm_index);
	
	rhs_matrix(problem_domain_mesh.mesh_cells(e, i), problem_domain_mesh.mesh_cells(f, j)) += ErichsenQuadRuleFlat(@LaplaceSLPKernel3D, 1, test_function_space{i}, ansatz_function_space{j}, shape_function_space, shape_function_space, e, f, panel_distance_matrix, neighboring_type_matrix, problem_domain_mesh.mesh_cells, problem_domain_mesh.mesh_nodes, -problem_domain_mesh.cell_normal_vectors(e, :), -problem_domain_mesh.cell_normal_vectors(f, :), problem_domain_mesh.cell_surface_metrics(e), problem_domain_mesh.cell_surface_metrics(f), problem_domain_mesh.min_cell_range, sobolev_function_space_order, test_function_space_order, galerkin_estimate_norm_index);
      endfor
    endfor
    
    progress_counter += 1;
    fprintf(stderr(), cstrcat("Percentage: %%", num2str(floor(progress_counter / progress_total_steps * 100)), "\r"));
  endfor
endfor

## Calculate RHS vector.
rhs_vector = rhs_matrix * neumann_bc;

## Solve the BEM problem for \f$u\f$, the potential distribution on the domain
## surface.
u = stiffness_matrix \ rhs_vector;

## Calculate the L2 error.
err = norm(u - u_star) / norm(u_star)

## Plot results.
NewFigure;
hold on;
plot(u, "ro-");
plot(u_star, "b+-");
legend({"u (BEM)", "u (analytical)"});
hold off;
