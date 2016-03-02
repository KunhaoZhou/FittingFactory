////////////////////////////////////////
// This file contains fit_polynomial class
////////////////////////////////////////
#ifndef fit_polynomial_h_
#define fit_polynomial_h_
#include "fit_function.h"
using namespace std;

static double polynomial(vnl_vector<double> params, double xi) {
  double y = 0;
  for (unsigned i = 0; i < params.size(); i++)
    y += params[i] * pow(xi, i);
  return y;
}

// contains scatter matrix
class fit_polynomial_least_squares_fn : public fit_least_squares_function {

public:
	fit_polynomial_least_squares_fn(vnl_vector<double>& xi, vnl_vector<double>& y): fit_least_squares_function(xi, y, 4){}
	~fit_polynomial_least_squares_fn(){}
	virtual void f(vnl_vector<double> const& params, vnl_vector<double>& residuals)	{
    residuals.set_size(y_.size());
		for (unsigned i = 0; i < residuals.size(); ++i)
      residuals[i] = y_[i] - polynomial(params, xi_[i]);
	}

	//for polynomial, 1 x xx xxx
	virtual void gradf(vnl_vector<double> const& params, vnl_matrix<double>& jacobian) {
		unsigned n = params.size();
		unsigned m = y_.size();
		jacobian.set_size(m, n);
		for (unsigned r = 0; r<m; ++r) {
			jacobian[r][0] = -1;
			jacobian[r][1] = -xi_[r];
			jacobian[r][2] = -xi_[r] * xi_[r];
			jacobian[r][3] = -xi_[r] * xi_[r] * xi_[r];
		}
	}
};

//
class fit_polynomial :public fit_function
{
public:  
  fit_polynomial() {
    xi_ = 0, y_ = 0;
    lsqf_ = new fit_polynomial_least_squares_fn(xi_, y_);
  }
  ~fit_polynomial(){};
  virtual string name(){ return "polynomial"; }
  virtual string type(){ return "linear"; }

  virtual void initParams(vnl_vector<double>& params) {
    params.set_size(4);
    params[0] = 1;
    params[1] = 1;
    params[2] = 1;
    params[3] = 1;
  }

	double g(double xi) { 
    if (!ifnormal_)
      return polynomial(params_, xi);
		double xnorm = (xi - xaverage_) / xstd_;
    return polynomial(params_, xnorm) * ystd_ + yaverage_;
	}
};

#endif // fit_polynomial_h_
