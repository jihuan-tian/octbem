simple_corners = GenerateSimplexCorners(3);
tetra_tproduct(@(x, y, z) 1, 5, simple_corners(:, 1), simple_corners(:, 2), simple_corners(:, 3))

simple_corners(2, 1) = 2.0;
simple_corners(3, 2) = 3.0;
simple_corners(4, 3) = 4.0;
tetra_tproduct(@(x, y, z) 1, 5, simple_corners(:, 1), simple_corners(:, 2), simple_corners(:, 3))
