////////////////////////////////////////
// This file contains fit_lm_solver class
////////////////////////////////////////
#ifndef fit_lm_solver_h_
#define fit_lm_solver_h_
#include "fit_solver.h"
using namespace std;

class fit_lm_solver :public fit_non_linear_solver {
public:
  fit_lm_solver(){}
  ~fit_lm_solver(){}
  virtual string name(){ return "levenberg_marquardt"; }
  bool solve(fit_function* fn) {
    fn_ = fn;
    fn_->initParams(params_);
    vnl_levenberg_marquardt levmarq(*fn->lsqf());
    levmarq.set_max_function_evals(12);  //intentionally call retry function
    levmarq.minimize(params_);
    fn->importParams(params_);
    levmarq.diagnose_outcome();
    return true;
  }

  bool retry(unsigned n_iter, double xtol, double ftol) {
    fn_->initParams(params_);
    vnl_levenberg_marquardt levmarq(*fn_->lsqf());
    levmarq.set_verbose(true);
    levmarq.set_max_function_evals(n_iter);
    levmarq.set_x_tolerance(xtol);
    levmarq.set_f_tolerance(ftol);
    levmarq.minimize(params_);
    fn_->importParams(params_);
    levmarq.diagnose_outcome();
    return true;
  };

private:
  fit_function* fn_;
};


#endif // fit_lm_solver_h_

