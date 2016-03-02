////////////////////////////////////////
// This file contains fit_factory class
////////////////////////////////////////

#ifndef fit_factory_h_
#define fit_factory_h_
#pragma comment(lib, "v3p_netlib.lib")
#include <unordered_map>
#include "fit_polynomial.h"
#include "fit_gaussian.h"
#include "fit_exponential.h"
#include "fit_svd_solver.h"
#include "fit_lm_solver.h"
#include "fit_amoeba_solver.h"
#include "fit_radial_basis_function.h"

class fit_factory {
public:
  fit_factory(){}
  ~fit_factory(){}
  fit_function* fn();
  fit_solver* fs();
  bool register_solver(fit_solver* fs);
  bool register_function(fit_function* fn);
  bool fit(vnl_vector<double> xi, vnl_vector<double> y, string fns, string fss, vnl_vector<double>& params);

private:
  fit_solver* fs_;
  fit_function* fn_;
  unordered_map<string, fit_function*> register_function_;
  unordered_map<string, fit_solver*> register_solver_;
};

#endif // fit_factory_h_
