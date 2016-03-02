////////////////////////////////////////
// This file contains svd solver test
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;
void test_svd_solver() {
	fit_factory ffac;
	ffac.register_solver(new fit_svd_solver);
	ffac.register_function(new fit_polynomial);
	unsigned m = 1000;
  vnl_vector<double> xi(m), y(m), params(4), fit_params(4);
  params[0] = -4; params[1] = -3; params[2] = 2; params[3] = 1;
	for (unsigned i = 0; i < m; i++) { 
    xi[i] = i; 
    y[i] = polynomial(params, xi[i]);
  }
  ffac.fit(xi, y, "polynomial", "svd", fit_params);
	cout << "y(xi=3) = " << y[5] << "<---original value" << endl;
	cout << "g(xi=3) = " << ffac.fn()->g(5) << "<---polynomial solved by svd" << endl << endl;
	system("pause");
}

TESTMAIN(test_svd_solver);