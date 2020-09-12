function pullback_kernel_values = BEMKernelPullback(xi, eta, nx_functor, ny_functor, kernel_function, kx_shape_functions_for_geometry, ky_shape_functions_for_geometry, kx_cell_node_coord_list, ky_cell_node_coord_list)
  ## BEMKernelPullback - Pullback a BEM kernel to the reference cell.
  ## @param xi A list of area coordinates on the reference cell \f$\hat{K}_x\f$.
  ## @param eta A list of area coordinates on the reference cell \f$\hat{K}_y\f$.
  ## @param nx_functor The functor depending on area coordinates for
  ## calculating the cell normal vectors for \f$K_x\f$.
  ## @param ny_functor The functor depending on area coordinates for
  ## calculating the cell normal vectors for \f$K_x\f$.
  ## @param kx_shape_functions_for_geometry A list of shape functions on the reference cell
  ## \f$\hat{K}_x\f$ for describing the real cell geometry.
  ## @param ky_shape_functions_for_geometry A list of shape functions on the reference cell
  ## \f$\hat{K}_y\f$ for describing the real cell geometry.
  ## @param kx_cell_node_coord_list A list of node global coordinates for the
  ## cell \f$K_x\f$.
  ## @param ky_cell_node_coord_list A list of node global coordinates for the
  ## cell \f$K_y\f$.

  number_of_area_coords = size(xi, 1);
  if (number_of_area_coords != size(eta, 1))
    error("The total numbers of the source and field point area coordinates should be the same!");
  else
    ## Transform the area coordinates \f$\xi\f$ in \f$\hat{K}_x\f$ and
    ## \f$\eta\f$ in \f$\hat{K}_y\f$ to global coordinates \f$x\f$ in \f$K_x\f$
    ## and \f$y\f$ in \f$K_y\f$.
    x = AreaToGlobalCoords(xi, kx_shape_functions_for_geometry, kx_cell_node_coord_list);
    y = AreaToGlobalCoords(eta, ky_shape_functions_for_geometry, ky_cell_node_coord_list);
    ## Calculate the normal vectors located at \f$\xi\f$.
    nx = nx_functor(xi);
    ## Calculate the normal vectors located at \f$\eta\f$;
    ny = ny_functor(eta);

    ## Substitute the global coordinates \f$x\f$ and \f$y\f$ into the
    ## kernel function. N.B. The normal vector data may not be used or
    ## fully used by the kernel function.
    pullback_kernel_values = kernel_function(x, y, nx, ny);
  endif
endfunction
