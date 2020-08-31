ConfigGraphicsToolkit;

enable_drawing = false;
enable_printing = false;

h = 1;
figure(h);
clf(h);

hold on;

## Plot the normalized Legendre polynomials from order 0 to 6.
x = -1:1e-2:1;
list_of_orders = 0:6;
for m = list_of_orders
  plot(x, LegendrePN(x, m), "linewidth", fig_line_width, "color", clist(m + 1, :));
endfor
grid on;
xlabel("X");
ylabel("Y");
title("Legendre polynomials");
legend(arrayfun(@(x) cstrcat("n=", num2str(x)), list_of_orders, "UniformOutput", false), "location", "EastOutside");

hold off;

if (enable_printing)
  PrintGCF("legendre-polynomials");
endif

## Generate the 1D Gauss-Legendre quadrature abscissas and weights.
for m = 1:10
  fprintf(stderr(), "=== Number of quadrature points: %d\n", m);
  [qpts, qwts] = GaussLegendreRule(m)
endfor

## Generate the 2D Gauss-Legendre quadrature abscissas and weights in the 2D square.
ref_rect = [-1, -1, 2, 2];
[qpts_square, qwts_square] = GaussLegendreRuleOnSquare([10, 10]);

h = 2;
figure(h);
clf(h);

hold on;

Plot2DRect(ref_rect, "k-");
Plot2DPointList(qpts_square, "k.", 20);

hold off;

grid on;
axis equal;
xlabel("X");
ylabel("Y");
title("Gauss-Legendre quadrature points");

if (enable_printing)
  PrintGCF("gauss-legendre-quadrature-points-in-square");
endif

## Generate the 3D Gauss-Legendre quadrature abscissas and weights in the 3D cube.
[qpts_cube, qwts_cube] = GaussLegendreRuleOnCube([10, 6, 3]);

h = 3;
figure(h);
clf(h);

hold on;

Plot3DPointList(qpts_cube, "ro");

hold off;

grid on;
axis equal;
xlabel("X");
ylabel("Y");
title("Gauss-Legendre quadrature points");

if (enable_printing)
  PrintGCF("gauss-legendre-quadrature-points-in-cube");
endif

## Generate Gauss-Jacobi related rules.
figure; Plot3DPointList(GaussJacobiRuleOnCube(0.5, 0.5, [10, 6, 3]), "ro");
figure; Plot2DPointList(GaussJacobiRuleOnSquare(0.5, 0.5, [10, 6]), "ro"); axis equal; xlim([-1, 1]); ylim([-1, 1]);
figure; Plot2DPointList(GaussJacobiLegendreRuleOnSquare(0.5, 0.5, [6, 6]), "ro"); axis equal; xlim([-1, 1]); ylim([-1, 1]);
figure; Plot2DPointList(GaussJacobiLegendreRuleOnSquare(0.5, 0.5, [10, 6]), "ro"); axis equal; xlim([-1, 1]); ylim([-1, 1]);
figure; Plot2DPointList(GaussLegendreJacobiRuleOnSquare(0.5, 0.5, [6, 6]), "ro"); axis equal; xlim([-1, 1]); ylim([-1, 1]);
figure; Plot2DPointList(GaussLegendreJacobiRuleOnSquare(0.5, 0.5, [10, 6]), "ro"); axis equal; xlim([-1, 1]); ylim([-1, 1]);

## Apply the quad rule.
ApplyQuadRule(@sin, qpts, qwts)
