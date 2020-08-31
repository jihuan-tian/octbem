function ret = inner_prod(v1, v2)
  ## \ingroup RadarCalib
  ## inner_prod - Calculate the inner product of two vectors.
  ## @param v1 the first input vector.
  ## @param v2 the second input vector.
  ## @param ret Inner product of the two vectors.
  
  ret = sum(v1 .* v2);
endfunction
