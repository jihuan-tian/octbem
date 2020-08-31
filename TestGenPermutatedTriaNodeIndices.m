tria_node_indices_matrix = GenTriaNodeIndicesMatrix(1:15)

for m = 1:3
  GenPermutatedTriaNodeIndicesForward(m, tria_node_indices_matrix)
endfor

for m = 1:3
  GenPermutatedTriaNodeIndicesBackward(m, tria_node_indices_matrix)
endfor
