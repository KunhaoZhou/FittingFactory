////////////////////////////////////////
// This file contains fit_amoeba_solver and fit_cost_function class
////////////////////////////////////////
#ifndef fit_amoeba_solver_h_
#define fit_amoeba_solver_h_
#include "fit_solver.h"
#include <vnl/algo/vnl_amoeba.h>
#include <vnl/vnl_cost_function.h>
using namespace std;

// The cost function (somewhat like a functor). 
class fit_cost_function : public vnl_cost_function {
public:
  fit_cost_function(fit_least_squares_function* f) :
    vnl_cost_function(f->get_number_of_unknowns()), f_(f), fx_(vnl_vector<double>(f->get_number_of_residuals())){}
	~fit_cost_function(){}
	double f(vnl_vector<double> const& x) {
		f_->f(x, fx_);
		double c = fx_.squared_magnitude();
		return c;
	}
protected:  
	fit_least_squares_function* f_;
	vnl_vector<double> fx_;
};


class fit_amoeba_solver :public fit_non_linear_solver {
public:
	virtual string name() { return "amoeba"; }
  virtual bool solve(fit_function* fn) {
		fn_ = fn;
    fn_->initParams(params_);
    fit_cost_function fcf(fn->lsqf());
		vnl_amoeba amoeba(fcf);
		//amoeba.set_x_tolerance(1e-14);
		//amoeba.set_f_tolerance(1e-14);
		amoeba.set_max_iterations(500);
    amoeba.minimize(params_);
    fn->importParams(params_);
		return true;
	}

  virtual bool retry(unsigned n_iter, double xtol, double ftol) {
    fn_->initParams(params_);
    fit_cost_function fcf(fn_->lsqf());
		vnl_amoeba amoeba(fcf);
    //amoeba.verbose = true;		
    //amoeba.set_relative_diameter(1.0);
		amoeba.set_x_tolerance(xtol);
		amoeba.set_f_tolerance(ftol);
		amoeba.set_max_iterations(n_iter);
    amoeba.minimize(params_);
    fn_->importParams(params_);
		return true;
	}
private:
	fit_function* fn_;
};

#endif // fit_amoeba_solver_h_

