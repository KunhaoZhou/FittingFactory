////////////////////////////////////////
// This file contains fit_radial_basis_function class
////////////////////////////////////////
#ifndef fit_radial_basis_function_h_
#define fit_radial_basis_function_h_
#include "fit_function.h"
using namespace std;

static double radial_basis(vnl_vector<double> w, double b, vnl_vector<double> chi, double x) {
	unsigned n = w.size();
	double s = 0;
	for (unsigned i = 0; i < n-1; i++) {
		double d = (x - chi[i])*(x - chi[i]);
		s = s + w[i] * exp(-b*d);
	}
	s = s + w[n - 1];  //add another parameter
	return s;
}

//contains f(x,residual) gradf(x,jacobian)
class fit_radial_least_squares_fn : public fit_least_squares_function {
public:
  fit_radial_least_squares_fn(vnl_vector<double>& xi, vnl_vector<double>& y) :
    fit_least_squares_function(xi, y, 4){ chi_.set_size(3); chi_[0] = 1; chi_[1] = 2; chi_[2] = 3; b_ = 1.5; }
	~fit_radial_least_squares_fn(){}

	void f(vnl_vector<double> const& params, vnl_vector<double>& residuals) {
    residuals.set_size(y_.size());
    for (unsigned i = 0; i < residuals.size(); ++i) {
			residuals[i] = y_[i] - radial_basis(params, b_, chi_, xi_[i]);
		}
	}
	void gradf(vnl_vector<double> const& params, vnl_matrix<double>& jacobian) {
		unsigned n = params.size();
		unsigned m = y_.size();
		jacobian.set_size(m, n);
		for (unsigned r = 0; r<m; ++r) {
			jacobian[r][0] = -exp(-b_*(xi_[r] - chi_[0])*(xi_[r] - chi_[0]));
			jacobian[r][1] = -exp(-b_*(xi_[r] - chi_[1])*(xi_[r] - chi_[1]));
			jacobian[r][2] = -exp(-b_*(xi_[r] - chi_[2])*(xi_[r] - chi_[2]));
			jacobian[r][3] = -1;
		}
	}

private:
	double b_;
	vnl_vector<double> chi_;
};

class fit_radial_basis_function :public fit_function {
public:
	fit_radial_basis_function() {
    xi_ = 0, y_ = 0;
    lsqf_ = new fit_radial_least_squares_fn(xi_, y_);
    b_ = 1.5;
    chi_.set_size(3);
    chi_[0] = 1;
    chi_[1] = 2;
    chi_[2] = 3;
  }
	~fit_radial_basis_function(){}
  virtual string name(){ return "radial_basis"; }

  virtual string type(){ return"non-linear"; }

  virtual void initParams(vnl_vector<double>& params) {
    params.set_size(4);
    params[0] = 1;
    params[1] = 1;
    params[2] = 1;
    params[3] = 1;
  }

	virtual double g(double x) {
		if (ifnormal_)
      return (ystd_*radial_basis(params_, b_, chi_, (x - xaverage_) / xstd_) + yaverage_);
		else
      return radial_basis(params_, b_, chi_, x);
	}
private:
	double b_;
	vnl_vector<double> chi_;
};

#endif // fit_radial_basis_function_h_
