ConfigGraphicsToolkit;

tria_order = 10;
tria_area_coords = AreaCoordsOnTriangle(tria_order);
std_tria_node_indices = GenTriaNodeIndicesMatrix(1:NumberOfNodesOnTriangle(tria_order));

## Permutation by starting from the 1st corner node in the forward direction.
NewFigure; Plot2DPointList(tria_area_coords(GenPermutatedTriaNodeIndicesForward(1, std_tria_node_indices), 1:2), "ro-");
## Permutation by starting from the 2nd corner node in the forward direction.
NewFigure; Plot2DPointList(tria_area_coords(GenPermutatedTriaNodeIndicesForward(2, std_tria_node_indices), 1:2), "ro-");
## Permutation by starting from the 3rd corner node in the forward direction.
NewFigure; Plot2DPointList(tria_area_coords(GenPermutatedTriaNodeIndicesForward(3, std_tria_node_indices), 1:2), "ro-");

## Permutation by starting from the 1st corner node in the backward direction.
NewFigure; Plot2DPointList(tria_area_coords(GenPermutatedTriaNodeIndicesBackward(1, std_tria_node_indices), 1:2), "ro-");
## Permutation by starting from the 2nd corner node in the backward direction.
NewFigure; Plot2DPointList(tria_area_coords(GenPermutatedTriaNodeIndicesBackward(2, std_tria_node_indices), 1:2), "ro-");
## Permutation by starting from the 3rd corner node in the backward direction.
NewFigure; Plot2DPointList(tria_area_coords(GenPermutatedTriaNodeIndicesBackward(3, std_tria_node_indices), 1:2), "ro-");
