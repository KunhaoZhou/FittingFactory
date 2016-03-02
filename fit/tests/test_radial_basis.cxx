////////////////////////////////////////
// This file contains test for radial basis function 
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;

static void test_radial_basis() {
	fit_factory ffac;
	ffac.register_solver(new fit_amoeba_solver);
	ffac.register_solver(new fit_lm_solver);
	ffac.register_function(new fit_radial_basis_function);

	unsigned m = 4;
  vnl_vector<double> xi(m), y(m), fit_params(m);
	vnl_vector<double> w(3+1), chi(3);
	chi[0] = 1; chi[1] = 2; chi[2] = 3;
	w[0] = 1; w[1] = 2; w[2] = 3; w[3] = 4;
  for (unsigned i = 0; i < m; i++) {
    xi[i] = i;
    y[i] = radial_basis(w, 1.5, chi, xi[i]); 
  }

	//levenberg_marquardt solver
  ffac.fit(xi, y, "radial_basis", "levenberg_marquardt", fit_params);
	cout << "y(xi=2) = " << y[2] << "<---original value" << endl;
	cout << "g(xi=2) = " << ffac.fn()->g(2) << "<---radial basis function solved by levenberg_marquardt" << endl;
	system("pause");

	//amoeba solver
  ffac.fit(xi, y, "radial_basis", "amoeba", fit_params);
	cout << "y(xi=2) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=2) = " << ffac.fn()->g(3) << "<---radial basis function solved by amoeba" << endl;
	system("pause");
}

TESTMAIN(test_radial_basis);