////////////////////////////////////////
// This file contains test for amoeba solver
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;

void test_amoeba_solver() {
	fit_factory ffac;
	ffac.register_solver(new fit_amoeba_solver);
	ffac.register_function(new fit_gaussian);
	ffac.register_function(new fit_exponential);
	ffac.register_function(new fit_polynomial);
	ffac.register_function(new fit_radial_basis_function);

	//gaussian test
	unsigned m = 9;
  vnl_vector<double> xi(m), y(m), fit_params;
	for (unsigned i = 0; i < m; i++) {
    xi[i] = i;		
    y[i] = gauss(10, 2, 5, 10, xi[i]);
  }
  ffac.fit(xi, y, "gaussian", "amoeba", fit_params);
	cout << "y(xi=3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=3) = " << ffac.fn()->g(3) << "<---gaussian solved by amoeba" << endl << endl;
	system("pause");


	//exponential test
	for (unsigned i = 0; i < m; i++) {
    xi[i] = i*0.1;
		y[i] = 2 * exp(0.25*xi[i]) + 2; 
  }
  ffac.fit(xi, y, "exponential", "amoeba", fit_params);
	cout << "y(xi=0.3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=0.3) = " << ffac.fn()->g(0.3) << "<---exponential solved by amoeba" << endl;
	system("pause");


	//polynomial test
	vnl_vector<double>	xx(4);
	xx[0] = -4; xx[1] = -3; xx[2] = 2; xx[3] = 1;
	for (unsigned i = 0; i < m; i++) {
    xi[i] = i; 
    y[i] = polynomial(xx, xi[i]);
  }
  ffac.fit(xi, y, "polynomial", "amoeba", fit_params);
	cout << "y(xi=1) = " << y[1] << "<---original value" << endl;
	cout << "g(xi=1) = " << ffac.fn()->g(1) << "<---polynomial solved by amoeba" << endl << endl;
	system("pause");

	//radial basis
	vnl_vector<double> w(3), chi(3);// xxi(6), yy(6);
	chi[0] = 1; chi[1] = 2; chi[2] = 3;
	w[0] = 10; w[1] = 5; w[2] = 2;
	for (unsigned i = 0; i < m; i++) { 
    xi[i] = i;
    y[i] = radial_basis(w, 1.5, chi, xi[i]);
  }

  fit_params.set_size(3);
  ffac.fit(xi, y, "radial_basis", "amoeba", fit_params);
	cout << "y(xi=2) = " << y[2] << "<---original value" << endl;
	cout << "g(xi=2) = " << ffac.fn()->g(2) << "<---radial basis solved by amoeba" << endl;
	system("pause");
}

TESTMAIN(test_amoeba_solver);