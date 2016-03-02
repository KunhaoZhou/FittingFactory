////////////////////////////////////////
// This file contains fit_function class
////////////////////////////////////////

#ifndef fit_function_h_
#define fit_function_h_
#include <iostream>
#include <vcl_vector.h>
#include <vnl/vnl_vector.h>
#include "fit_least_squares_function.h"

using namespace std;
static double sd(const vnl_vector<double> x); // standard deviation
static void getnormal(vnl_vector<double>& x); // normalization

class fit_function {
public:	
	virtual string name() = 0;
	virtual string type() = 0; // linear or non-linear
	virtual double g(double x) = 0; // optimized model
  virtual void initParams(vnl_vector<double>& x) = 0;
  virtual void importData(vnl_vector<double> xi, vnl_vector<double> y);
  void importParams(vnl_vector<double> params) { params_ = params; } // import result from solver
  fit_least_squares_function* lsqf();
  vnl_vector<double> xi();
  vnl_vector<double> y();
  vnl_vector<double> params();
  void clear();
  void normalize();
protected:
  fit_least_squares_function* lsqf_;
	vnl_vector<double> xi_;
  vnl_vector<double> y_;
  vnl_vector<double> params_; // model parameters
	double xaverage_;
	double yaverage_;
	double xstd_;
  double ystd_;
  bool ifnormal_;
};

#endif // fit_function_h_
