/////////////////////////////////////////
// This file contains test for exponential 
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;

static void test_exponential() {
	fit_factory ffac;
	ffac.register_solver(new fit_amoeba_solver);
	ffac.register_solver(new fit_lm_solver);
	ffac.register_function(new fit_exponential);

	unsigned m = 100;
  vnl_vector<double> xi(m), y(m), fit_params;
	//exponential test
	for (unsigned i = 0; i < m; i++) {
		xi[i] = 0.1*i;
		y[i] = 2 * exp(0.25*xi[i]) + 2;
	}
  ffac.fit(xi, y, "exponential", "amoeba", fit_params);
	cout << "y(xi=0.3) = " << y[3]<<"<---original value" << endl;
	cout << "g(xi=0.3) = " << ffac.fn()->g(0.3) <<"<---exponential solved by amoeba"<< endl;
	//ffac.fs()->retry(1000000, 1e-15, 1e-15);
	system("pause");

  ffac.fit(xi, y, "exponential", "levenberg_marquardt", fit_params);
	cout << "y(xi=0.2) = " << y[2] << "<---original value" << endl;
	cout << "g(xi=0.2) = " << ffac.fn()->g(0.2) << "<---exponential solved by levenberg_marquardt" << endl;
	system("pause");
}

TESTMAIN(test_exponential);