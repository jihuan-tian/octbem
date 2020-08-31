function is_aligned = IsVectorListAligned(vector_list)
  ## IsVectorListAligned - Determine if a list of vectors in 3D space are
  ## aligned to the same direction.
  ## @param vector_list A list of vectors.
  ## @param is_aligned Whether the vectors are aligned.

  [number_of_vectors, space_dim] = size(vector_list);
  
  for m = 1:(number_of_vectors - 1)
    for n = (m + 1):number_of_vectors
      if (abs(inner_prod(NormalizeVect(vector_list(m, :)), NormalizeVect(vector_list(n, :))) - 1) > 1e-10)
	is_aligned = false;
	
	return;
      endif
    endfor
  endfor
  
  is_aligned = true;
endfunction
