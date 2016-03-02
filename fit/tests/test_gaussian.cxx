////////////////////////////////////////
// This file contains test for gaussian 
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;

static void test_gaussian()
{
	fit_factory ffac;
	ffac.register_solver(new fit_amoeba_solver);
	ffac.register_solver(new fit_lm_solver);
	ffac.register_function(new fit_gaussian);

	//gaussian test
	unsigned m = 9;
  vnl_vector<double> xi(m), y(m), fit_params;
	for (unsigned i = 0; i < m; i++)	{ xi[i] = i;		y[i] = gauss(10, 2, 5, 10, xi[i]); }
  ffac.fit(xi, y, "gaussian", "levenberg_marquardt", fit_params);
	cout << "y(xi=3) = " << y[3] << "<---original value"<<endl;
	cout << "g(xi=3) = " << ffac.fn()->g(3) << "<---gaussian solved by levenberg_marquardt" << endl << endl;
	system("pause");

  ffac.fit(xi, y, "gaussian", "amoeba", fit_params);
	cout << "y(xi=4) = " << y[4] << "<---original value" << endl;
	cout << "g(xi=4) = " << ffac.fn()->g(4) << "<---gaussian solved by amoeba" << endl << endl;
	system("pause");
}
TESTMAIN(test_gaussian);