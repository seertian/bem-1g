#include "Matrial.h"

Matrial::Matrial() {}
Matrial::Matrial(const int g) {
	group_num = g;
	upper_energy_boundary.resize(g+1);
	d.resize(g);
	siga.resize(g);
	nusigf.resize(g);
	sigt.resize(g);
	sign2n.resize(g);
	chi.resize(g);
	scattering_matrix.resize(g,g);
}
Matrial::~Matrial() {}
///////// put vector
void Matrial::put_d(VecDoub &xs) {
	d = xs;
}
void Matrial::put_siga(VecDoub &xs) {
	siga = xs;
}
void Matrial::put_nusigf(VecDoub &xs) {
	nusigf = xs;
}
void Matrial::put_sigt(VecDoub &xs ) {
	sigt = xs;
}
void Matrial::put_sign2n(VecDoub &xs) {
	sign2n = xs;
}
void Matrial::put_chi(VecDoub &xs) {
	chi = xs;
}
void Matrial::put_scattering(MatDoub &xs) {
	scattering_matrix = xs;
}
//////////////////// put signal value
void Matrial::put_d(int g, double  &xs) {
	d[g] = xs;
}
void Matrial::put_siga(int g, double  &xs) {
	siga[g] = xs;
}
void Matrial::put_nusigf(int g, double  &xs) {
	nusigf[g] = xs;
}
void Matrial::put_sigt(int g, double  &xs ) {
	sigt[g] = xs;
}
void Matrial::put_sign2n(int g, double  &xs) {
	sign2n[g] = xs;
}
void Matrial::put_chi(int g, double  &xs) {
	chi[g] = xs;
}
void Matrial::put_scattering(int i, int j, double  &xs) {
	scattering_matrix[i][j] = xs;
}
/////////////////////////// get vector

VecDoub Matrial::get_d() {
	return d;
}
VecDoub Matrial::get_siga() {
	return siga;
}
VecDoub Matrial::get_nusigf() {
	return nusigf;
}
VecDoub Matrial::get_sigt( ) {
	return sigt;
}
VecDoub Matrial::get_sign2n() {
	return sign2n;
}
VecDoub Matrial::get_chi() {
	return chi;
}
MatDoub Matrial::get_scattering() {
	return scattering_matrix;
}
////////////////////
double Matrial::get_d(int g) {
	return d[g];
}
double Matrial::get_siga(int g) {
	return siga[g];
}
double Matrial::get_nusigf(int g) {
	return nusigf[g];
}
double Matrial::get_sigt(int g) {
	return sigt[g];
}
double Matrial::get_sign2n(int g) {
	return sign2n[g];
}
double Matrial::get_chi(int g) {
	return chi[g];
}
double Matrial::get_scattering(int i, int j) {
	return scattering_matrix[i][j];
}
///////////////////////////
void Matrial::print_xs() {

}
SimpleMatrial::SimpleMatrial() {}

SimpleMatrial::SimpleMatrial(double xs1, double xs2, double xs3) {
	d = xs1;
	siga = xs2;
	nusigf = xs3;
}

SimpleMatrial::~SimpleMatrial() {}

double SimpleMatrial::get_d() {
	return this->d;
}

double SimpleMatrial::get_siga() {
	return this->siga;
}

double SimpleMatrial::get_nusigf() {
	return this->nusigf;
}
void SimpleMatrial::put_d(double xs) {
	d = xs;
}
void SimpleMatrial::put_siga(double xs) {
	siga = xs;
}
void SimpleMatrial::put_nusigf(double xs) {
	nusigf = xs;
}
double SimpleMatrial::diffusion() {
	return this->d;
}

double SimpleMatrial::absorption() {
	return this->siga;
}

double SimpleMatrial::production() {
	return this->nusigf;
}

