function ret = average(v, varargin)
  ## average - Calculate the average of a vector or matrix.
  ##
  ## If the input argument v is a matrix, by default, the average operation is
  ## performed for each column. If an additional string argument "r" is
  ## provided, the average operation is carried out for rows. Any other string
  ## will be interpreted as average in column.
  
  [m, n] = size(v);

  ## If the input data is a vector.
  if (m == 1) || (n == 1)
    ret = sum(v) / length(v);
  else
    if (length(varargin) == 1)
      switch (varargin{1})
	case 'r'
	  ## Average the data in rows.
	  ret = zeros(m, 1);

	  for i = 1:m
	    ret(i) = sum(v(i,:)) / n;
	  endfor
	otherwise
	  ## Average the data in columns.
	  ret = zeros(1, n);

	  for i = 1:n
	    ret(i) = sum(v(:,i)) / m;
	  endfor
      endswitch
    else
      ## Average the data in columns.
      ret = zeros(1, n);

      for i = 1:n
	ret(i) = sum(v(:,i)) / m;
      endfor
    endif
  endif
endfunction
