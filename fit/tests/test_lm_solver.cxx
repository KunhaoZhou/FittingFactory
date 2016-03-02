////////////////////////////////////////
// This file contains test for levenberg_marquardt solver
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;

void test_lm_solver() {
	fit_factory ffac;
	ffac.register_solver(new fit_lm_solver);
	ffac.register_function(new fit_gaussian);
	ffac.register_function(new fit_exponential);
	ffac.register_function(new fit_polynomial);
	ffac.register_function(new fit_radial_basis_function);

	//gaussian test
	unsigned m = 100;
  vnl_vector<double> xi(m), y(m), fit_params;
	for (unsigned i = 0; i < m; i++)	{ xi[i] = i;		y[i] = gauss(10, 2, 5, 0, xi[i]); }
  ffac.fit(xi, y, "gaussian", "levenberg_marquardt", fit_params);
	cout << "y(xi=4) = " << y[4] << "<---original value" << endl;
	cout << "g(xi=4) = " << ffac.fn()->g(4) << "<---gaussian solved by levenberg_marquardt" << endl << endl;
	system("pause");


	//exponential test
	for (unsigned i = 0; i < m; i++)
	{
		xi[i] = i*0.1;
		y[i] = 2 * exp(0.25*xi[i]) + 2;
	}
  ffac.fit(xi, y, "exponential", "levenberg_marquardt", fit_params);
	cout << "y(xi=0.3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=0.3) = " << ffac.fn()->g(0.3) << "<---exponential solved by levenberg_marquardt" << endl;
	system("pause");


	//polynomial test
	vnl_vector<double>	params(4);
  params[0] = -4; params[1] = -3; params[2] = 2; params[3] = 1;
	for (unsigned i = 0; i < m; i++) {
    xi[i] = i;
    y[i] = polynomial(params, xi[i]);
  }
  ffac.fit(xi, y, "polynomial", "levenberg_marquardt", fit_params);
	cout << "y(xi=4) = " << y[4] << "<---original value" << endl;
	cout << "g(xi=4) = " << ffac.fn()->g(4) << "<---polynomial solved by levenberg_marquardt" << endl << endl;
	system("pause");

	//radial basis
	vnl_vector<double> w(3), chi(3);
	chi[0] = 1; chi[1] = 2; chi[2] = 3;
	w[0] = 10; w[1] = 5; w[2] = 2;
	for (unsigned i = 0; i < m; i++) { 
    xi[i] = i; 
    y[i] = radial_basis(w, 1.5, chi, xi[i]); 
  }

  fit_params.set_size(3);
  ffac.fit(xi, y, "radial_basis", "levenberg_marquardt", fit_params);
	cout << "y(xi=3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=3) = " << ffac.fn()->g(3) << "<---radial basis solved by levenberg_marquardt" << endl;
	system("pause");
}
TESTMAIN(test_lm_solver);