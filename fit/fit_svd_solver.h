////////////////////////////////////////
// This file contains fit_svd_solver class
////////////////////////////////////////
#ifndef fit_svd_solver_h_
#define fit_svd_solver_h_
#include "fit_solver.h"
#include "fit_function.h"
#include <vnl/algo/vnl_svd.h>
using namespace std;

//product of vectors
static vnl_vector<double> dotp(vnl_vector<double> v1, vnl_vector<double> v2) {
	vnl_vector<double> v3(v1.size());
	for (unsigned i = 0; i < v1.size() && v1.size() == v2.size(); i++)
		v3[i] = v1[i] * v2[i];
	return v3;
}

//expectation matrix
static vnl_matrix<double> expmatrix(vnl_vector<double> xi, vnl_vector<double> y) {
	vnl_matrix<double> EM(5, 5);
	vnl_vector<double> xx, xxx, xxxx, xxxxx, xxxxxx, xy, xxy, xxxy, yy;
	xx = dotp(xi, xi);
	xxx = dotp(xx, xi);
	xxxx = dotp(xx, xx);
	xxxxx = dotp(xxx, xx);
	xxxxxx = dotp(xxx, xxx);
	xy = dotp(xi, y);
	xxy = dotp(xx, y);
	xxxy = dotp(xxx, y);
	yy = dotp(y, y);
	unsigned m = y.size();
	EM[0][0] = xxxxxx.mean(); EM[0][1] = xxxxx.mean(); EM[0][2] = xxxx.mean(); EM[0][3] = xxxy.mean(); EM[0][4] = xxx.mean();
	EM[1][0] = xxxxx.mean(); EM[1][1] = xxxx.mean(); EM[1][2] = xxx.mean(); EM[1][3] = xxy.mean(); EM[1][4] = xx.mean();
	EM[2][0] = xxxx.mean(); EM[2][1] = xxx.mean(); EM[2][2] = xx.mean(); EM[2][3] = xy.mean(); EM[2][4] = xi.mean();
	EM[3][0] = xxxy.mean(); EM[3][1] = xxy.mean(); EM[3][2] = xy.mean(); EM[3][3] = yy.mean(); EM[3][4] = y.mean();
	EM[4][0] = xxx.mean(); EM[4][1] = xx.mean(); EM[4][2] = xi.mean(); EM[4][3] = y.mean(); EM[4][4] = 1;
	return EM;
}


class fit_svd_solver :public fit_linear_solver {
public:
	fit_svd_solver(){}
	~fit_svd_solver(){}
	virtual string name(){ return "svd"; }

	bool solve(fit_function* fn) {
		fn_ = fn;
    fn_->initParams(params_);
    xi_ = fn->lsqf()->xi();
		y_ = fn->lsqf()->y();  //get normalized data
		M_ = expmatrix(xi_, y_); //get expectation matrix
		vnl_svd<double> svd(M_); //use svd to solve M
		svd.zero_out_absolute();
		V_ = svd.V();
    params_ = V_.get_column(4);
    double c = -params_[3];
    for (unsigned i = 0; i < params_.size(); i++)
      params_[i] = params_[i] / c;
		vnl_vector<double> params(4);
    params[0] = params_[4];
    params[1] = params_[2];
    params[2] = params_[1];
    params[3] = params_[0];      //in the correct order
    params_ = params;
    fn->importParams(params);
		return true;
	}

	bool retry(unsigned n_iter, double xtol, double ftol)	{ return true; }

private:
	fit_function* fn_; 
	vnl_matrix<double> M_;
	vnl_matrix<double> V_;
	vnl_vector<double> xi_;
	vnl_vector<double> y_;
};

#endif // fit_svd_solver_h_

