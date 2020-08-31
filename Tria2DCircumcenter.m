function circumcenter = Tria2DCircumcenter(tria_corners)
  ## Tria2DCircumcenter - Calculate the circumcenter for a triangle in 2D space.
  ## @param tria_corners The 2D coordinates for triangle corners, stored as a 3*2 matrix.
  ## @param circumcenter The calculated 2D coordinate for the circumcenter.
  
  circumcenter(1) = -(tria_corners(2, 2)*(tria_corners(3, 2)^2-tria_corners(1, 2)^2+tria_corners(3, 1)^2-tria_corners(1, 1)^2)+tria_corners(1, 2)*(-tria_corners(3, 2)^2-tria_corners(3, 1)^2+tria_corners(2, 1)^2)+tria_corners(1, 2)^2*tria_corners(3, 2)-tria_corners(2, 1)^2*tria_corners(3, 2)+tria_corners(1, 1)^2*tria_corners(3, 2)+tria_corners(2, 2)^2*(tria_corners(1, 2)-tria_corners(3, 2)))/(2*tria_corners(2, 1)*tria_corners(3, 2)-2*tria_corners(1, 1)*tria_corners(3, 2)+(2*tria_corners(1, 1)-2*tria_corners(3, 1))*tria_corners(2, 2)+(2*tria_corners(3, 1)-2*tria_corners(2, 1))*tria_corners(1, 2));
  circumcenter(2) = (tria_corners(2, 1)*(tria_corners(3, 2)^2+tria_corners(3, 1)^2-tria_corners(1, 1)^2)+tria_corners(1, 1)*(-tria_corners(3, 2)^2-tria_corners(3, 1)^2)+(tria_corners(1, 1)-tria_corners(3, 1))*tria_corners(2, 2)^2+(tria_corners(3, 1)-tria_corners(2, 1))*tria_corners(1, 2)^2+tria_corners(1, 1)^2*tria_corners(3, 1)+tria_corners(2, 1)^2*(tria_corners(1, 1)-tria_corners(3, 1)))/(2*tria_corners(2, 1)*tria_corners(3, 2)-2*tria_corners(1, 1)*tria_corners(3, 2)+(2*tria_corners(1, 1)-2*tria_corners(3, 1))*tria_corners(2, 2)+(2*tria_corners(3, 1)-2*tria_corners(2, 1))*tria_corners(1, 2));
endfunction
