#ifndef _LINALG_H_
#define _LINALG_H_
#include "DataStructure.h"
// this the a simple linalg libaray with function operation

MatDoub lu_decomposition(MatDoub &a) {
	double TINY  =1.0e-40;
	
	double _temp;
	int i, j, k, i_max;
	double temp, big_tmp;
	int n = a.nrows();
	MatDoub lu(n, n), a_ref(n, n);
	lu = a, a_ref = a;
	VecDoub vv(n), index(n);
	double d = 1.0;
	for (i = 0; i < n; i++) {
		big_tmp =0.0;
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
				_temp = lu[i_max][j];
				lu[i_max][j] = lu[k][j];
				lu[k][j] = _temp;
			}
			d = -d;
			vv[i_max] = vv[k];
		}
		index[k] = i_max;
		if (lu[k][k] == 0.0) lu[k][k] = TINY;
		for (i = k + 1; i < n; i++) {
			_temp = lu[i][k] /= lu[k][k];
			for (j = k + 1; j < n; j++)
				lu[i][j] -= _temp * lu[k][j];
		}
	}
	MatDoub lu_with_index(n + 1, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			lu_with_index[i + 1][j] = lu[i][j];
		}
		lu_with_index[0][i] = index[i];
	}
	return lu_with_index;
}

VecDoub solve_with_lu(MatDoub &a, VecDoub &b) {
	int n = a.nrows();
	int ii = 0, ip;
	VecDoub index(n);
	VecDoub x(n);
	double sum;
	MatDoub lu(n, n), lu_with_index;
	lu_with_index = lu_decomposition(a);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			lu[i][j] = lu_with_index[i + 1][j];
		}
		index[i] = lu_with_index[0][i];
	}
	//VecDoub x = b;

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

MatComplex lu_decomposition(MatComplex &a) {
	Complex TINY  (1.0e-40,1.0e-40);
	
	Complex _temp;
	int i, j, k, i_max;
	double temp, big_tmp;
	int n = a.nrows();
	MatComplex lu(n, n), a_ref(n, n);
	lu = a, a_ref = a;
	VecDoub vv(n), index(n);
	double d = 1.0;
	for (i = 0; i < n; i++) {
		big_tmp =0.0;
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
				_temp = lu[i_max][j];
				lu[i_max][j] = lu[k][j];
				lu[k][j] = _temp;
			}
			d = -d;
			vv[i_max] = vv[k];
		}
		index[k] = i_max;
		if (lu[k][k] == 0.0) lu[k][k] = TINY;
		for (i = k + 1; i < n; i++) {
			_temp = lu[i][k] /= lu[k][k];
			for (j = k + 1; j < n; j++)
				lu[i][j] -= _temp * lu[k][j];
		}
	}
	MatComplex lu_with_index(n + 1, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			lu_with_index[i + 1][j] = lu[i][j];
		}
		lu_with_index[0][i] = index[i];
	}
	return lu_with_index;
}

VecComplex solve_with_lu(MatComplex &a, VecComplex &b) {
	int n = a.nrows();
	int ii = 0, ip;
	VecDoub index(n);
	VecComplex x(n);
	Complex sum;
	MatComplex lu(n, n), lu_with_index;
	lu_with_index = lu_decomposition(a);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			lu[i][j] = lu_with_index[i + 1][j];
		}
		index[i] = lu_with_index[0][i].real();
		// this is special for complex number 
	}
	//VecDoub x = b;

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
// For QR
template <class T>
T norm(Vector<T> &x) {
	T sum;
	for (int i = 0; i < x.size(); ++i)
	{
		sum+=x[i]*x[i];
	}
	return sqrt(sum);
}
// For QR
template <class T>
T scalar_product(Vector<T> &x, Vector<T> &y) {
    T sum;
    int len=x.size();
    for (int i = 0; i < len; i++) {
        sum += x[i]*y[i];
    }
    return sum;
}

// find the projection of x on y
template <class T>
void projection(Vector<T> &x, Vector<T> &y, int len) {
	//int len =x.size();
	T sp = scalar_product(x,y);
	T norm_y2 = norm(y)*norm(y);
	for (int i = 0; i < len; ++i)
	{
		x[i]=x[i]-sp/norm_y2*y[i];
	}
}
// template <class T>
// Matrix<T> qr_decomposition(Matrix<T> &a) {
// 	int m = a.nrows();
// 	int n = a.ncols();

// 	Matrix<T> v(n, m);
// 	for (int i = 0; i < m; i++) {
// 		for (int j = 0; j < n; j++) {
// 			v[j][i] = a[i][j];
// 		}
// 	}
// 	for (int i = 1; i < n; i++) {
// 		for (int j = 0; j < i; j++) {
// 			projection(v[i], v[j], m);
// 		}
// 	}

// 	Matrix<T> q(m, n);
// 	q = v.transpose();	
// 	q = q.normalize_by_col();
	
// 	Matrix<T> q_transposed;
// 	q_transposed = q.transpose();

// 	Matrix<T> r(n, n);
// 	for (int i = 0; i < n; i++) {
//         for (int j = 0;j < n;j++) {
//             r[i][j] = 0.0;
//             for (int k = 0; k < m; k++) {
//                 r[i][j] += q_transposed[i][k] * a[k][j];
//             }
//         }
//     }

//     Matrix<T> qr(m+n, n);
//     for (int i = 0; i < m; ++i) {
//     	for (int j = 0; j < n; ++j) {
//     		qr[i][j] = q[i][j];
//     	}
//     }
//     for (int i = 0; i < n; ++i) {
//     	for (int j = 0; j < n; ++j) {
//     		qr[i+m][j] = r[i][j];
//     	}
//     }
//     Matrix<T> zero_matrix(m+n, n);
// 	qr = zero_matrix - qr;
//     return qr;
// }
#endif /* _LINALG_H_ */