////////////////////////////////////////
// This file contains fit_solver class
////////////////////////////////////////
#ifndef fit_solver_h_
#define fit_solver_h_
#include "fit_function.h"
#include <iostream>

class fit_solver {
public:
	fit_solver(){}
	~fit_solver(){}
	virtual bool solve(fit_function* fn) = 0;
  virtual string name() = 0;
  virtual string type() = 0;	
  virtual bool retry(unsigned n_iter, double xtol, double ftol) = 0;
  vnl_vector<double> params() { return params_; }
protected:
  vnl_vector<double> params_;
};


class fit_non_linear_solver :public fit_solver {
public:
	virtual string type(){ return "non-linear"; }
};

class fit_linear_solver :public fit_solver {
public:
  virtual string type(){ return "linear"; }
};

#endif // fit_solver_h_
