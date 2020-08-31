function v = NormalizeVect(varargin)
  ## NormalizeVect - Normalize a vector.
  ##
  ## varargin{1}: the vector.
  ##
  ## varargin{2}: optional, the norm number p.

  eps_effective_line_seg = 1e-10;
  
  if (nargin == 2)
    p = varargin{2};
  else
    p = 2;
  endif

  vector_norm = norm(varargin{1}, p);
  if (vector_norm < eps_effective_line_seg)
    v = zeros(size(varargin{1}));
  else
    v = varargin{1} / vector_norm;
  endif
endfunction
