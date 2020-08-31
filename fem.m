
ConfigGraphicsToolkit;

## Nodes in the distorted triangle.
tri_node_list_distorted = [-4.01406, 9.94323, 4.5;
                           -3.65318, 6.16233, 3.8;
                           -2.29768, 2.279, 3;
                           -0.282742, -0.981539, 3.8;
                           2.8038, -4.13172, 4.5;
                           -6.69391, 6.71186, 3.8;
                           -5.41167, 2.64535, 3;
                           -3.14029, -1.5677, 3;
                           -0.575822, -4.82824, 3.8;
                           -8.41576, 3.48796, 3;
                           -6.40083, -0.725094, 2;
                           -3.87299, -4.82824, 3;
                           -9.58809, 0.00761196, 2;
                           -7.13353, -4.24207, 1;
                           -10.0638, -3.29188, 1];

## Generate and plot a triangle from area coordinate.
xi = 1:-0.01:0;
eta = 1 - xi;
zet = zeros(size(xi));
tri_edge1 = GenDistortedTria(xi, eta, zet);

eta = xi;
xi = zet;
zet = 1 - eta;
tri_edge2 = GenDistortedTria(xi, eta, zet);

zet = eta;
eta = xi;
xi = 1 - zet;
tri_edge3 = GenDistortedTria(xi, eta, zet);

figure;
Plot3DPointList([tri_edge1, tri_edge2, tri_edge3]', "r-");
hold on;
Plot3DPointList(tri_node_list_distorted, "bo");
xlabel("x");
ylabel("y");
zlabel("z");
PrintGCF("distorted-tria-from-local-coord");
