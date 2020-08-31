format long;

K1 = [3, 5, 0;
      6, 7, 0;
      8, 3, 0];
K2 = [4, 1, 0;
      9, -2, 0;
      7, -8, 0];

## The answer is: 3.342516087186934
distance = DistOfCoplanarConvexPolygons(K1, K2)
