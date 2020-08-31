# octbem
Boundary element method (BEM) demonstration using GNU Octave

Boundary element method (BEM) is an effective algorithm compared to finite element method (FEM) for resolving physical fields in open domain and/or complex models with geometric details, especially those having large dimensional scale difference. Typical application areas of BEM include electromagnetics, acoustics, crack propagation etc.

The basic idea of BEM is to construct the solution of a partial differential equation (PDE), like the 2nd order Laplace equation, by using a representation formula derived from the Green's identity. By approaching this representation formula to the domain boundary with some presumption on potential continuity, boundary integral equation can be obtained.

This project aims to demonstrate key technologies involved in Galerkin BEM. It is mainly developed in GNU Octave for easy-of-use. Meanwhile, those functions involving heavy computations are also gradually being replaced by compiled C++ object modules for improving efficiency. 

Existing functionalities in this project include handling of numerical quadrature for various scenarios, namely, regular, near singular, singular and hypersingular integrals. The method is verified by solving a Laplace problem with Neumann boundary condition.

For the theoretical background of the adopted method, please refer to: Erichsen, Stefan, and Stefan A. Sauter. 1998. “Efficient Automatic Quadrature in 3-d Galerkin BEM.” Computer Methods in Applied Mechanics and Engineering, Papers presented at the Seventh Conference on Numerical Methods and Computational Mechanics in Science and Engineering, 157 (3): 215–24.



