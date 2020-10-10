function lagrange_basis_functors = LagrangeBasisOn3DQuad(order)
  ## LagrangeBasisOn3DQuad - Generate a list of Lagrange basis
  ## functions constructed on the reference quadrangle depending on the
  ## two independent area coordinates. At present, the maximum
  ## supported Lagrange basis function order is 3. These basis
  ## functions are implemented as function handles.
  ## @param order The basis function order.
  ## @param lagrange_basis_functors The list of generated Lagrange basis functions.
  
  switch(order)
    case 1
      lagrange_basis_functors{1} = @(chi) (1-chi(:,2)).*(1-chi(:,1));
      lagrange_basis_functors{2} = @(chi) chi(:,2).*(1-chi(:,1));
      lagrange_basis_functors{3} = @(chi) (1-chi(:,2)).*chi(:,1);
      lagrange_basis_functors{4} = @(chi) chi(:,2).*chi(:,1);
    case 2
      lagrange_basis_functors{1} = @(chi) 4.*(chi(:,2)-1).*(chi(:,2)-1/2).*(chi(:,1)-1).*(chi(:,1)-1/2);
      lagrange_basis_functors{2} = @(chi) -8.*(chi(:,2)-1).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-1/2);
      lagrange_basis_functors{3} = @(chi) 4.*(chi(:,2)-1/2).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-1/2);
      lagrange_basis_functors{4} = @(chi) -8.*(chi(:,2)-1).*(chi(:,2)-1/2).*(chi(:,1)-1).*chi(:,1);
      lagrange_basis_functors{5} = @(chi) 16.*(chi(:,2)-1).*chi(:,2).*(chi(:,1)-1).*chi(:,1);
      lagrange_basis_functors{6} = @(chi) -8.*(chi(:,2)-1/2).*chi(:,2).*(chi(:,1)-1).*chi(:,1);
      lagrange_basis_functors{7} = @(chi) 4.*(chi(:,2)-1).*(chi(:,2)-1/2).*(chi(:,1)-1/2).*chi(:,1);
      lagrange_basis_functors{8} = @(chi) -8.*(chi(:,2)-1).*chi(:,2).*(chi(:,1)-1/2).*chi(:,1);
      lagrange_basis_functors{9} = @(chi) 4.*(chi(:,2)-1/2).*chi(:,2).*(chi(:,1)-1/2).*chi(:,1);
    case 3
      lagrange_basis_functors{1} = @(chi) (81.*(chi(:,2)-1).*(chi(:,2)-2/3).*(chi(:,2)-1/3).*(chi(:,1)-1).*(chi(:,1)-2/3).*(chi(:,1)-1/3))/4;
      lagrange_basis_functors{2} = @(chi) -(243.*(chi(:,2)-1).*(chi(:,2)-2/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-2/3).*(chi(:,1)-1/3))/4;
      lagrange_basis_functors{3} = @(chi) (243.*(chi(:,2)-1).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-2/3).*(chi(:,1)-1/3))/4;
      lagrange_basis_functors{4} = @(chi) -(81.*(chi(:,2)-2/3).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-2/3).*(chi(:,1)-1/3))/4;
      lagrange_basis_functors{5} = @(chi) -(243.*(chi(:,2)-1).*(chi(:,2)-2/3).*(chi(:,2)-1/3).*(chi(:,1)-1).*(chi(:,1)-2/3).*chi(:,1))/4;
      lagrange_basis_functors{6} = @(chi) (729.*(chi(:,2)-1).*(chi(:,2)-2/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-2/3).*chi(:,1))/4;
      lagrange_basis_functors{7} = @(chi) -(729.*(chi(:,2)-1).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-2/3).*chi(:,1))/4;
      lagrange_basis_functors{8} = @(chi) (243.*(chi(:,2)-2/3).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-2/3).*chi(:,1))/4;
      lagrange_basis_functors{9} = @(chi) (243.*(chi(:,2)-1).*(chi(:,2)-2/3).*(chi(:,2)-1/3).*(chi(:,1)-1).*(chi(:,1)-1/3).*chi(:,1))/4;
      lagrange_basis_functors{10} = @(chi) -(729.*(chi(:,2)-1).*(chi(:,2)-2/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-1/3).*chi(:,1))/4;
      lagrange_basis_functors{11} = @(chi) (729.*(chi(:,2)-1).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-1/3).*chi(:,1))/4;
      lagrange_basis_functors{12} = @(chi) -(243.*(chi(:,2)-2/3).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-1).*(chi(:,1)-1/3).*chi(:,1))/4;
      lagrange_basis_functors{13} = @(chi) -(81.*(chi(:,2)-1).*(chi(:,2)-2/3).*(chi(:,2)-1/3).*(chi(:,1)-2/3).*(chi(:,1)-1/3).*chi(:,1))/4;
      lagrange_basis_functors{14} = @(chi) (243.*(chi(:,2)-1).*(chi(:,2)-2/3).*chi(:,2).*(chi(:,1)-2/3).*(chi(:,1)-1/3).*chi(:,1))/4;
      lagrange_basis_functors{15} = @(chi) -(243.*(chi(:,2)-1).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-2/3).*(chi(:,1)-1/3).*chi(:,1))/4;
      lagrange_basis_functors{16} = @(chi) (81.*(chi(:,2)-2/3).*(chi(:,2)-1/3).*chi(:,2).*(chi(:,1)-2/3).*(chi(:,1)-1/3).*chi(:,1))/4;
    otherwise
      lagrange_basis_functors = {};
      error("At present, the maximum supported Lagrange basis function order is 3.");
  endswitch
endfunction
