////////////////////////////////////////
// This file contains fit_least_squares_function class
////////////////////////////////////////
#ifndef fit_least_squares_function_h
#define fit_least_squares_function_h
#include <vnl/algo/vnl_levenberg_marquardt.h>
#include <vnl/vnl_least_squares_function.h>
class fit_least_squares_function :public vnl_least_squares_function {
public:
  fit_least_squares_function(vnl_vector<double> const& xi, vnl_vector<double> const& y, unsigned int n) :
    vnl_least_squares_function(n, xi.size(), vnl_least_squares_function::use_gradient),
    xi_(xi), y_(y) {}
  virtual ~fit_least_squares_function(){};
  virtual void f(vnl_vector<double> const& x, vnl_vector<double>& fx) = 0;
  virtual void gradf(vnl_vector<double> const& x, vnl_matrix<double>& jacobian) = 0;
  vnl_vector<double> xi(){ return xi_; }
  vnl_vector<double> y(){ return y_; }
  void importData(vnl_vector<double>& xi, vnl_vector<double>& y){ xi_ = xi;	y_ = y; this->n_ = y_.size(); }
  void clear(){ xi_ = 0;	y_ = 0; }

protected:
  vnl_vector<double> xi_;
  vnl_vector<double> y_;
};

#endif
