#ifndef __MATRIAL_H_
#define __MATRIAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream> // read matrix 
#include <iomanip> // print format control
#include <ccomplex>
#include "DataStructure.h" //Vector and Matrix Data Stucture

using namespace std;

class Matrial {
private:
	int group_num;
	int n = group_num + 1;
	VecDoub upper_energy_boundary;
	/////////////////////////////
	VecDoub d;
	VecDoub siga;
	VecDoub nusigf;

	VecDoub sigt;
	VecDoub sign2n;
	VecDoub chi;
	////////////////////////////
	MatDoub scattering_matrix;

public:
	Matrial();
	Matrial(int);
	//Matrial(double, double, double);
	//Matrix(char *);
	~Matrial();
	void put_d(VecDoub &);
	void put_siga(VecDoub &);
	void put_nusigf(VecDoub &);
	void put_sigt(VecDoub &);
	void put_sign2n(VecDoub &);
	void put_chi(VecDoub &);
	void put_scattering(MatDoub &);

	void put_d(int, double &);
	void put_siga(int, double &);
	void put_nusigf(int, double &);
	void put_sigt(int, double &);
	void put_sign2n(int, double &);
	void put_chi(int, double &);
	void put_scattering(int, int, double &);

	VecDoub get_d();
	VecDoub get_siga();
	VecDoub get_nusigf();
	VecDoub get_sigt();
	VecDoub get_sign2n();
	VecDoub get_chi();
	MatDoub get_scattering();

	double get_d(int);
	double get_siga(int);
	double get_nusigf(int);
	double get_sigt(int);
	double get_sign2n(int);
	double get_chi(int);
	double get_scattering(int, int);

	void print_xs();
};
class SimpleMatrial {
private:
	//int group_num;
	//int n=group+1;
	//VecDoub upper_energy_boundary(n,0.0);
	double d;
	double siga;
	double nusigf;

	double sigt;
	double sign2n;
	double chi;
	MatDoub scattering_matrix;

public:
	SimpleMatrial();
	SimpleMatrial(double, double, double);
	//SimpleMatrial(char *);
	~SimpleMatrial();

	double get_d();
	double get_siga();
	double get_nusigf();
	double diffusion();
	double absorption();
	double production();
	void put_d(double);
	void put_siga(double);
	void put_nusigf(double);
	// void put_scattering(MatDoub);

};
#endif /* defined(__MATRIAL_H_) */