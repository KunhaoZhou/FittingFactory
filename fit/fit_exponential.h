////////////////////////////////////////
// This file contains fit_exponential class
////////////////////////////////////////
#ifndef fit_exponential_h_
#define fit_exponential_h_
#include "fit_function.h"
using namespace std;

//contains f(x,residual) gradf(x,jacobian)
class fit_exponential_least_squares_fn : public fit_least_squares_function {
public:
	fit_exponential_least_squares_fn(vnl_vector<double>& xi, vnl_vector<double>& y) :fit_least_squares_function(xi, y, 3){}
	~fit_exponential_least_squares_fn(){}

  virtual void f(vnl_vector<double> const& params, vnl_vector<double>& residuals) {
		residuals = y_;
		double alpha = params[0];
		double b = params[1];
		double c = params[2];
		for (unsigned i = 0; i < residuals.size(); ++i) {
			residuals[i] = y_[i] - alpha*exp(b*xi_[i]) - c;
		}
	}

  virtual void gradf(vnl_vector<double> const& params, vnl_matrix<double>& jacobian) {
		unsigned n = params.size();
		unsigned m = y_.size();
		jacobian.set_size(m, n);
		double alpha = params[0];
		double b = params[1];
		double c = params[2];
		for (unsigned r = 0; r<m; ++r) {
			jacobian[r][0] = -exp(b*xi_[r]);
			jacobian[r][1] = -(alpha*xi_[r])* exp(b*xi_[r]);
			jacobian[r][2] = -1;
		}
	}
};

//
class fit_exponential:public fit_function {
public:
	fit_exponential() {
    xi_ = 0, y_ = 0;
    lsqf_ = new fit_exponential_least_squares_fn(xi_, y_);
  }
	~fit_exponential(){}
	virtual string name(){ return "exponential"; }
  virtual string type(){ return"non-linear"; }
  virtual void initParams(vnl_vector<double>& params) {
    params.set_size(3);
    params[0] = 1;
    params[1] = 1;
    params[2] = 0;
  }

	virtual double g(double x) {
		if (ifnormal_)
      return ystd_*(params_[0] * exp(params_[1] * (x - xaverage_) / xstd_) + params_[2]) + yaverage_;
		else
      return (params_[0] * exp(params_[1] * x) + params_[2]);
	}
};

#endif // fit_exponential_h_

