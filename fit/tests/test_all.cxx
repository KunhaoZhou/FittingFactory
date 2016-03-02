////////////////////////////////////////
// This file contains test all problem solutions
////////////////////////////////////////
#include <testlib/testlib_test.h>
#include "fit/fit_factory.h"
#include <iostream>
using namespace std;

void test_all() {   
  //registration
	fit_factory ffac;
	ffac.register_solver(new fit_amoeba_solver);
	ffac.register_solver(new fit_lm_solver);
	ffac.register_solver(new fit_svd_solver);
	ffac.register_function(new fit_gaussian);
	ffac.register_function(new fit_exponential);
	ffac.register_function(new fit_polynomial);
	ffac.register_function(new fit_radial_basis_function);
	unsigned m = 100;
  vnl_vector<double> xi(m), y(m), fit_params;

	//svd + polynomial test
  vnl_vector<double>	params(4);
  params[0] = -4; params[1] = -3; params[2] = 2; params[3] = 1;
	for (unsigned i = 0; i < m; i++) { 
    xi[i] = i; 
    y[i] = polynomial(params, xi[i]);
  }
  ffac.fit(xi, y, "polynomial", "svd", fit_params);
	cout << "y(xi=3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=3) = " << ffac.fn()->g(3) << "<---polynomial solved by svd" << endl << endl;
	system("pause");


	//amoeba + gaussian test

	for (unsigned i = 0; i < m; i++) { 
    xi[i] = i*2;		
    y[i] = gauss(10, 2, 5, 10, xi[i]); 
  }
  ffac.fit(xi, y, "gaussian", "amoeba", fit_params);
	cout << "y(xi=6) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=6) = " << ffac.fn()->g(6) << "<---gaussian solved by amoeba" << endl << endl;
	system("pause");


	//amoeba + exponential test
	for (unsigned i = 0; i < m; i++) {
    xi[i] = i*0.1;		
    y[i] = 2 * exp(0.25*xi[i]) + 2;	
  }
  ffac.fit(xi, y, "exponential", "amoeba", fit_params);
	cout << "y(xi=0.3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=0.3) = " << ffac.fn()->g(0.3) << "<---exponential solved by amoeba" << endl;
	system("pause");


	//amoeba + polynomial test
	for (unsigned i = 0; i < m; i++) {
    xi[i] = i; 
    y[i] = polynomial(params, xi[i]);
  }
  ffac.fit(xi, y, "polynomial", "amoeba", fit_params);
	cout << "y(xi=1) = " << y[1] << "<---original value" << endl;
	cout << "g(xi=1) = " << ffac.fn()->g(1) << "<---polynomial solved by amoeba" << endl << endl;
	system("pause");

	//amoeba + radial basis
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


	//lm + gaussian test
  for (unsigned i = 0; i < m; i++) {
    xi[i] = i;		
    y[i] = gauss(10, 2, 5, 10, xi[i]);
  }
  ffac.fit(xi, y, "gaussian", "levenberg_marquardt", fit_params);
	cout << "y(xi=4) = " << y[4] << "<---original value" << endl;
	cout << "g(xi=4) = " << ffac.fn()->g(4) << "<---gaussian solved by levenberg_marquardt" << endl << endl;
	system("pause");


	//lm + exponential test
	for (unsigned i = 0; i < m; i++) {
    xi[i] = i*0.1;
    y[i] = 2 * exp(0.25*xi[i]) + 2;
  }
  ffac.fit(xi, y, "exponential", "levenberg_marquardt", fit_params);
	cout << "y(xi=0.3) = " << y[3] << "<---original value" << endl;
	cout << "g(xi=0.3) = " << ffac.fn()->g(0.3) << "<---exponential solved by levenberg_marquardt" << endl;
	system("pause");


	//lm + polynomial test
  params[0] = -4; params[1] = -3; params[2] = 2; params[3] = 1;
	for (unsigned i = 0; i < m; i++) { 
    xi[i] = i; 
    y[i] = polynomial(params, xi[i]);
  }
  ffac.fit(xi, y, "polynomial", "levenberg_marquardt", fit_params);
	cout << "y(xi=4) = " << y[4] << "<---original value" << endl;
	cout << "g(xi=4) = " << ffac.fn()->g(4) << "<---polynomial solved by levenberg_marquardt" << endl << endl;
	system("pause");

	//lm + radial basis
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

TESTMAIN(test_all);