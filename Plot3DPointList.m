function ret = Plot3DPointList(varargin)
  ## Plot3DPointList - Plot a list of 3D points.
  ##
  ## varargin: the first argument is a Nx3 matrix. Each row represents a 3D
  ## point and the columns represents XYZ components. The second optional
  ## argument specifies the line type and/or point type for the plotting.
  ##
  ## ret: the figure handle.
  
  ret = [];
  
  if (length(varargin) >= 2)
    pt_style = varargin{2};
  else
    pt_style = "ro";
  endif
  
  ret = plot3(varargin{1}(:,1), varargin{1}(:,2), varargin{1}(:,3), pt_style);
endfunction
