////////////////////////////////////////
// This file contains test for polynomial 
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;

static void test_polynomial() {
	fit_factory ffac;
	ffac.register_solver(new fit_amoeba_solver);
	ffac.register_solver(new fit_lm_solver);
	ffac.register_solver(new fit_svd_solver);
	ffac.register_function(new fit_polynomial);

	unsigned m = 10;
  vnl_vector<double> xi(m), y(m), fit_params;
	// polynomial test
	vnl_vector<double>	params(4);
  params[0] = -4; params[1] = -3; params[2] = 2; params[3] = 1;
  for (unsigned i = 0; i < m; i++) {
    xi[i] = i; 
    y[i] = polynomial(params, xi[i]);
  }

  ffac.fit(xi, y, "polynomial", "levenberg_marquardt", fit_params);
	cout << "y(xi=3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=3) = " << ffac.fn()->g(3) << "<---polynomial solved by levenberg_marquardt" << endl << endl;
	system("pause");

  ffac.fit(xi, y, "polynomial", "amoeba", fit_params);
	cout << "y(xi=3) = " << y[4] << "<---original value" << endl;
	cout << "g(xi=3) = " << ffac.fn()->g(4) << "<---polynomial solved by amoeba" << endl << endl;
	system("pause");


  ffac.fit(xi, y, "polynomial", "svd", fit_params);
	cout << "y(xi=3) = " << y[5] << "<---original value" << endl;
	cout << "g(xi=3) = " << ffac.fn()->g(5) << "<---polynomial solved by svd" << endl << endl;
	system("pause");

}

TESTMAIN(test_polynomial);