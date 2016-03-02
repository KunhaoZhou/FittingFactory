#include "fit_function.h"

static double sd(const vnl_vector<double> x) {
  double a = x.mean(), s = 0;
  for (unsigned i = 0; i < x.size(); i++) {
    s = s + (x[i] - a)*(x[i] - a);
  }
  return sqrt(s / x.size());
}

static void getnormal(vnl_vector<double>& x) {
  double a = x.mean();
  double s = sd(x);
  for (unsigned i = 0; i < x.size(); i++)
    x[i] = (x[i] - a) / s;
}

fit_least_squares_function* fit_function::lsqf() { 
  return lsqf_; 
}

vnl_vector<double> fit_function::xi() {
  return xi_; 
}

vnl_vector<double> fit_function::y() {
  return y_; 
}

vnl_vector<double> fit_function::params(){
  return params_; 
}

void fit_function::normalize() {
  vnl_vector<double> norxi, nory;
  norxi = xi_;
  nory = y_;
  xaverage_ = xi_.mean();
  yaverage_ = y_.mean();
  xstd_ = sd(xi_);
  ystd_ = sd(y_);
  getnormal(norxi);
  getnormal(nory);
  lsqf_->importData(norxi, nory);
  ifnormal_ = 1; //to decide which g(x) to adopt
}

void fit_function::clear() {
  ifnormal_ = 0;
  xi_ = 0;
  y_ = 0;
  params_ = 0;
  lsqf_->clear();
}

void fit_function::importData(vnl_vector<double> xi, vnl_vector<double> y) {
  xi_ = xi;
  y_ = y;
  lsqf_->importData(xi, y);
  //this->normalize();
}
