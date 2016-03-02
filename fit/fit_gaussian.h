////////////////////////////////////////
// This file contains fit_gaussian class
////////////////////////////////////////

#ifndef fit_gaussian_h_
#define fit_gaussian_h_
#include "fit_function.h"
using namespace std;

// gaussian model
static double gauss(double alpha, double mu, double sigma, double c, double xi) {
  double arg = (xi - mu) / sigma;
  return alpha*exp(-arg*arg / 2) + c;
}

// contains f(x,residual) gradf(x,jacobian)
class fit_gauss_least_squares_fn : public fit_least_squares_function {
public:
  fit_gauss_least_squares_fn(vnl_vector<double>& xi, vnl_vector<double>& y) :fit_least_squares_function(xi, y, 4){}
  ~fit_gauss_least_squares_fn(){}
  virtual void f(vnl_vector<double> const& gauss_params, vnl_vector<double>& residuals) {
    residuals.set_size(y_.size());
    double alpha = gauss_params[0];
    double mu = gauss_params[1];
    double sigma = gauss_params[2];
    double c = gauss_params[3];
    for (unsigned i = 0; i < residuals.size(); ++i)
      residuals[i] = y_[i] - gauss(alpha, mu, sigma, c, xi_[i]);
  }
  virtual void gradf(vnl_vector<double> const& gauss_params, vnl_matrix<double>& jacobian) {
    unsigned n = gauss_params.size();
    unsigned m = y_.size();
    jacobian.set_size(m, n);
    double alpha = gauss_params[0];
    double mu = gauss_params[1];
    double sigma = gauss_params[2];
    double sigma_2 = 1 / (sigma*sigma);
    double sigma_3 = sigma_2 / sigma;
    for (unsigned r = 0; r < m; ++r) {
      double d = xi_[r] - mu;
      jacobian[r][0] = -gauss(1.0, mu, sigma, 0, xi_[r]);
      jacobian[r][1] = -gauss(alpha, mu, sigma, 0, xi_[r])*d*sigma_2;
      jacobian[r][2] = -gauss(alpha, mu, sigma, 0, xi_[r])*d*d*sigma_3;
      jacobian[r][3] = -1;
    }
  }
};

//
class fit_gaussian :public fit_function {
public:
  fit_gaussian() {
    xi_ = 0;
    y_ = 0;
    lsqf_ = new fit_gauss_least_squares_fn(xi_, y_);
  }

  ~fit_gaussian(){};
  virtual string name(){ return "gaussian"; }
  virtual string type(){ return"non-linear"; }
  virtual double g(double x) {
    if (ifnormal_)
      return	yaverage_ + ystd_ * gauss(params_[0], params_[1], params_[2], params_[3], (x - xaverage_) / xstd_);
    return gauss(params_[0], params_[1], params_[2], params_[3], (x));
  }

  virtual void initParams(vnl_vector<double>& params) { //initial guass parameters
    params.set_size(4);
    double mean = 0, xsq = 0, ymax = -vnl_numeric_traits<double>::maxval;
    unsigned m = xi_.size();
    for (unsigned i = 0; i<m; ++i) {
      double temp = xi_[i];
      mean += temp;
      xsq += temp*temp;
      if (y_[i]>ymax) ymax = y_[i];
    }
    mean /= m;
    double var = xsq - m*mean*mean;
    params[0] = ymax;
    params[1] = mean;
    params[2] = sqrt(var / m);
    params[3] = 0;
  }
};

#endif // fit_gaussian_h_

