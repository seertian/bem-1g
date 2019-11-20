#include "LinAlg.h"

template <class T> //Complex or Doub/Double
Matrix<T> lu_decomposition(Matrix<T> &a) {
	double TINY = 1.0e-40;
	int i, j, k, i_max;
	double temp, big_tmp;
	int n = a.nrows();
	Matrix<T> lu(n, n), a_ref(n, n);
	lu = a, a_ref = a;
	Vector<T> vv(n), index(n);
	double d = 1.0;
	for (i = 0; i < n; i++) {
		big_tmp = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = abs(lu[i][j])) > big_tmp) big_tmp = temp;
		if (big_tmp == 0.0)
			cout << "Error, Singular matrix in LU decomposition" << endl;
		vv[i] = 1.0 / big_tmp;
	}
	for (k = 0; k < n; k++) {
		big_tmp = 0.0;
		for (i = k; i < n; i++) {
			temp = vv[i] * abs(lu[i][k]);
			if (temp > big_tmp) {
				big_tmp = temp;
				i_max = i;
			}
		}
		if (k != i_max) {
			for (j = 0; j < n; j++) {
				temp = lu[i_max][j];
				lu[i_max][j] = lu[k][j];
				lu[k][j] = temp;
			}
			d = -d;
			vv[i_max] = vv[k];
		}
		index[k] = i_max;
		if (lu[k][k] == 0.0) lu[k][k] = TINY;
		for (i = k + 1; i < n; i++) {
			temp = lu[i][k] /= lu[k][k];
			for (j = k + 1; j < n; j++)
				lu[i][j] -= temp * lu[k][j];
		}
	}
	Matrix<T> lu_with_index(n + 1, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			lu_with_index[i + 1][j] = lu[i][j];
		}
		lu_with_index[0][i] = index[i];
	}
	return lu_with_index;
}

template <class T> //Complex or Doub/Double
Vector<T> solve_with_lu(Matrix<T> &a, Vector<T> &b) {
	int n = a.nrows();
	int ii = 0, ip;
	vector<double> index(n);
	Vector<T> x(n);
	double sum;
	Matrix<T> lu(n, n), lu_with_index;
	lu_with_index = lu_decomposition(a);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			lu[i][j] = lu_with_index[i + 1][j];
		}
		index[i] = lu_with_index[0][i];
	}
	//Vector<T> x = b;

	for (int i = 0; i < n; i++) x[i] = b[i];

	for (int i = 0; i < n; i++) {
		ip = int(index[i]);
		sum = x[ip];
		x[ip] = x[i];
		if (ii != 0)
			for (int j = ii - 1; j < i; j++) sum -= lu[i][j] * x[j];
		else if (sum != 0.0) ii = i + 1;
		x[i] = sum;
	}
	for (int i = n - 1; i >= 0; i--) {
		sum = x[i];
		for (int j = i + 1; j < n; j++) sum -= lu[i][j] * x[j];
		x[i] = sum / lu[i][i];
	}
	return x;
}