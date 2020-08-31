#include <octave/oct.h>
#include <iterator>
#include <cmath>
#include <string>

#include "quadrule.hpp"

DEFUN_DLD (GaussJacobiRule, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn {} {} [qpts, qwts] = GaussJacobiRule (@var{order, alpha, beta})\n\
Calculate the Gauss-Jacobi quadrature points and weights on the standard domain [-1, 1].\n\
@end deftypefn\n\n\
@@param alpha Exponent in \f$(1-x)^{\alpha}\f$ multiplied to the integrand.\n\
@@param beta Exponent in \f$(1+x)^{\beta}\f$ multiplied to the integrand.\n\
@@param order The maximum polynomial order.\n\n\
@@param qpts The quadrature points.\n\
@@param qwts The quadrature weights.") {
  
  if (args.length() != 3) {
    print_usage();
  }

  //! The polynomial order to be generated.
  double alpha = args(0).double_value();
  double beta = args(1).double_value();
  unsigned int order = args(2).int_value();

  //! Allocate memory for the result matrix.
  Matrix qpts(1, order);
  Matrix qwts(1, order);

  jacobi_ek_compute(order, alpha, beta, (double *)qpts.fortran_vec(), (double *)qwts.fortran_vec());
  
  octave_value_list quad_rule;
  quad_rule(0) = octave_value(qpts);
  quad_rule(1) = octave_value(qwts);
  
  return quad_rule;
}
