#include <testlib/testlib_register.h>
DECLARE(test_all);
DECLARE(test_gaussian);
DECLARE(test_polynomial);
DECLARE(test_exponential);
DECLARE(test_radial_basis);
DECLARE(test_svd_solver);
DECLARE(test_lm_solver);
DECLARE(test_amoeba_solver);

void
register_tests() {
	REGISTER(test_all);
	REGISTER(test_gaussian);
  REGISTER(test_polynomial);
	REGISTER(test_exponential);
	REGISTER(test_radial_basis);
  REGISTER(test_svd_solver);	
  REGISTER(test_amoeba_solver);
	REGISTER(test_lm_solver);
};

DEFINE_MAIN;