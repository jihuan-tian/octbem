function neighboring_type_matrix = CalcQuadPanelNeighboringTypes(mesh_cells)
  ## CalcQuadPanelNeighboringTypes - Determine the neighboring type of each pair of panels.
  ## @param mesh_cell Cell data in the mesh.
  ## @param neighboring_type_matrix Neighboring type between each pair of panels.
  
  number_of_cells = size(mesh_cells, 1);
  neighboring_type_matrix = ones(number_of_cells, number_of_cells);

  for m = 1:number_of_cells
    for n = (m+1):number_of_cells
      neighboring_type_matrix(m, n) = GetQuadNeighboringType(m, n, mesh_cells);
      neighboring_type_matrix(n, m) = neighboring_type_matrix(m, n);
    endfor
  endfor
endfunction
