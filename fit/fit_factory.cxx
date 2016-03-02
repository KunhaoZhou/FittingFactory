#include "fit_factory.h"

fit_function* fit_factory::fn() { 
  return fn_; 
}

fit_solver* fit_factory::fs() {
  return fs_; 
}

bool fit_factory::register_solver(fit_solver* fs) {
  pair<string, fit_solver*> p(fs->name(), fs);
  register_solver_.insert(p);
  return true;
}

bool fit_factory::register_function(fit_function* fn) {
  pair<string, fit_function*> p(fn->name(), fn);
  register_function_.insert(p);
  return true;
}

bool fit_factory::fit(vnl_vector<double> xi, vnl_vector<double> y, string fns, string fss, vnl_vector<double>& params) {
  fs_ = register_solver_[fss]; // find solver according to name.
  fn_ = register_function_[fns]; // find function according to name.
  if ((!fn_) || (!fs_)) {
    cout << "Solution or Model not found" << endl;
    return false;
  }
  fn_->clear(); // clear data
  fn_->importData(xi, y);
  fs_->solve(fn_);   // solve
  params = fn_->params();   // extract result.

  vnl_vector<double> residuals;
  fn_->lsqf()->f(params, residuals);
  cout << endl << "Residuals sum = " << residuals.magnitude();
  if (residuals.magnitude() < (1e-10)) {
    cout << " <---Good enough. " << endl << endl;
    return true;
  }

  // need to retry. more iterations
  cout << "<--- Unsatisfactory first try. " << endl << endl << "Let's retry. " << endl;
  fs_->retry(100000, 1e-12, 1e-12);
  params = fn_->params();
  fn_->lsqf()->f(params, residuals);
  cout << endl << "Residuals sum = " << residuals.magnitude();
  if (residuals.magnitude() < (1e-5)) {
    cout << "  <---Good enough. " << endl << endl;
    return true;
  }

  // still unsatisfactory
  cout << "  <---Well...You might need to adopt another model or solver." << endl << endl;
  return false;
}