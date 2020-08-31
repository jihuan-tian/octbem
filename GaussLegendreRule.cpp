#include <octave/oct.h>
#include <iterator>
#include <cmath>
#include <string>

#include "quadrule.hpp"

DEFUN_DLD (GaussLegendreRule, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn {} {} GaussLegendreRule (@var{order})\n\
Calculate the Gauss-Legendre quadrature points and weights on the standard domain [-1, 1].\n\
@end deftypefn\n\n\
@@param order The maximum polynomial order.\n\n\
@@param qpts The quadrature points.\n\n\
@@param qwts The quadrature weights.") {

  if (args.length() != 1) {
    print_usage();
  }
  
  //! The polynomial order to be generated.
  unsigned int order = args(0).int_value();
  //! Allocate memory for the result matrix.
  Matrix qpts(order, 1);
  Matrix qwts(order, 1);

  legendre_dr_compute(order, (double *)qpts.fortran_vec(), (double *)qwts.fortran_vec());
  
  octave_value_list quad_rule;
  quad_rule(0) = octave_value(qpts);
  quad_rule(1) = octave_value(qwts);
  
  return quad_rule;
}
