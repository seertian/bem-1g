#include "DataStructure.h"

using namespace std; //including complex_literals

// this file contains the greens' function of Helmholtz Operator


// Complex gf1d(Complex &k, double &distance) {
// 	Complex gf = (exp(1i * k * abs(distance))) * 1i / (2.0 * k);
// 	return gf;
// }

double phi_star_core(double &distance, double &k)
{	
	double value = 0;
	value = -1/(2*k)*sin(k*abs(distance));
	return value;
}

double phi_star_diff_core(double &distance, double &k)
{	
	double value;
	value =  0.50*cos(k*distance);
	if (distance > 0) {
		value = -value;
	}
	return  value;
}

double phi_star_reflector(double &distance, double &k)
{	
	double value = 1/(2*k)*exp(-k*abs(distance));
	return value;
}

double phi_star_diff_reflector(double &distance, double &k)
{	
	double value = 0.50*exp(-k*distance);
	if (distance == 0) {		
		value = 0;
	}
	else if (distance > 0) {
		value = -value;
	}
	return value;
}

// for 3-d dimensional 

// Complex gf3d(Complex &k, double &distance) {
// 	Complex gf = (exp(-1i * k * abs(distance))) * 1i / (4.0*pi * abs(distance));
// 	return gf;
// }
