quad_node_indices_matrix = GenQuadNodeIndicesMatrix(1:16)

for m = 1:4
  GenPermutatedQuadNodeIndicesForward(m, quad_node_indices_matrix)
endfor

for m = 1:4
  GenPermutatedQuadNodeIndicesBackward(m, quad_node_indices_matrix)
endfor
