#include "LinAlg.h"
#include "EvaluationFunctionOneGroup.h"

template <class T> 
Matrix<T> jacobian_matrix(Vector<T> &x) {
	int n = x.size();
	Vector<T> fvec = EvaluationFunction(x);
	Matrix<T> fjac(n, n);
	Vector<T> xh = x;
	for (int j = 0; j < n; ++j) {
		double temp = xh[j];
		double h = EPS * abs(temp);
		if (h == 0.0) h = EPS;
		xh[j] = temp + h;
		h = xh[j] - temp;
		Vector<T> fxh = EvaluationFunction(xh);
		xh[j] = temp;
		for (int i = 0; i < n; ++i) {
			fjac[i][j] = (fxh[i] - fvec[i]) / h;
		}
		//return fjac;
	}
	return fjac;
}

template <class T> 
Vector<T> newton_raphson(Vector<T>  &ini_x) {
	const int ntrial = 1000;
	int n = ini_x.size();
	Vector<T>  x = ini_x;
	Vector<T>  p(n), fx(n);
	Matrix<T>  fjac(n, n);
	//
	for (int iter = 0; iter < ntrial; ++iter) {
		fx = EvaluationFunction(x);
		fjac = jacobian_matrix(x);
		double errf = 0.0;
		for (int i = 0; i < n; ++i) errf += abs(fx[i]);
		if (errf <= EPS) iter = ntrial + 1;

		for (int i = 0; i < n; ++i) p[i] = -fx[i];
		p = solve_with_lu(fjac, p);
		double errx = 0.0;
		for (int i = 0; i < n; ++i) {
			errx += abs(p[i]);
			x[i] += p[i];
		}
		if (errx <= EPS) iter = ntrial + 1;
		//old_x = x;
	}
	return x;
}

int main() {
	int n = 2;
	VecDoub x(n);
	x[0] = 1.4;
	x[1] = 0.1;

	MatDoub m(n, n);

	VecDoub sol(n);
	//sol[0] = exp(1);
	sol = newton_raphson(x);
	cout << "Root finding with two varibles(k_eff and neutron flux in core center) :" << endl;
	for (int i = 0; i < n; ++i) {
		cout << sol[i] << " " ;
	}
	cout << endl;
	VecDoub test(n);
	test = EvaluationFunction(sol);
	cout << "Error of Root finding: "<< endl;
	for (int i = 0; i < n; ++i) {
		cout << test[i] << " " ;
	}
	cout << endl;

}