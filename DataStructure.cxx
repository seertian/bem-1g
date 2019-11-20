#include "DataStructure.h"
// Vector definitions

template <class T>
Vector<T>::Vector() : nn(0), v(NULL) {}

template <class T>
Vector<T>::Vector(int n) : nn(n), v(n > 0 ? new T[n] : NULL) {}

template <class T>
Vector<T>::Vector(int n, const T& a) : nn(n), v(n > 0 ? new T[n] : NULL)
{
	for (int i = 0; i < n; i++) v[i] = a;
}

template <class T>
Vector<T>::Vector(int n, const T *a) : nn(n), v(n > 0 ? new T[n] : NULL)
{
	for (int i = 0; i < n; i++) v[i] = *a++;
}

template <class T>
Vector<T>::Vector(const Vector<T> &rhs) : nn(rhs.nn), v(nn > 0 ? new T[nn] : NULL)
{
	for (int i = 0; i < nn; i++) v[i] = rhs[i];
}

template <class T>
Vector<T> & Vector<T>::operator=(const Vector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != NULL) delete [] (v);
			nn = rhs.nn;
			v = nn > 0 ? new T[nn] : NULL;
		}
		for (int i = 0; i < nn; i++)
			v[i] = rhs[i];
	}
	return *this;
}

template <class T>
inline T & Vector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw ("Error! Vector subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline const T & Vector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw ("Error! Vector subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline int Vector<T>::size() const
{
	return nn;
}

template <class T>
void Vector<T>::resize(int newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

template <class T>
void Vector<T>::assign(int newn, const T& a)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i = 0; i < nn; i++) v[i] = a;
}

template <class T>
Vector<T>::~Vector()
{
	if (v != NULL) delete[] (v);
}

// end of Vector definitions



template <class T>
Matrix<T>::Matrix() : nn(0), mm(0), v(NULL) {}

template <class T>
Matrix<T>::Matrix(int n, int m) : nn(n), mm(m), v(n > 0 ? new T * [n] : NULL)
{
	int i, nel = m * n;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T &a) : nn(n), mm(m), v(n > 0 ? new T * [n] : NULL)
{
	int i, j, nel = m * n;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
	for (i = 0; i < n; i++) for (j = 0; j < m; j++) v[i][j] = a;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T *a) : nn(n), mm(m), v(n > 0 ? new T * [n] : NULL)
{
	int i, j, nel = m * n;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
	for (i = 0; i < n; i++) for (j = 0; j < m; j++) v[i][j] = *a++;
}

template <class T>
Matrix<T>::Matrix(const Matrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn > 0 ? new T * [nn] : NULL)
{
	int i, j, nel = mm * nn;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i, j, nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn = rhs.nn;
			mm = rhs.mm;
			v = nn > 0 ? new T*[nn] : NULL;
			nel = mm * nn;
			if (v) v[0] = nel > 0 ? new T[nel] : NULL;
			for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
		}
		for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
inline T* Matrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw ("Error! Matrix subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline const T* Matrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw ("Error! Matrix subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline int Matrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline int Matrix<T>::ncols() const
{
	return mm;
}

template <class T>
void Matrix<T>::resize(int newn, int newm)
{
	int i, nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn > 0 ? new T*[nn] : NULL;
		nel = mm * nn;
		if (v) v[0] = nel > 0 ? new T[nel] : NULL;
		for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
	}
}

template <class T>
void Matrix<T>::assign(int newn, int newm, const T& a)
{
	int i, j, nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn > 0 ? new T*[nn] : NULL;
		nel = mm * nn;
		if (v) v[0] = nel > 0 ? new T[nel] : NULL;
		for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
	}
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = a;
}
template <class T>
void Matrix<T>::print() {
	cout << setprecision(4);
	cout << setiosflags(ios::fixed);
	cout << setiosflags(ios::right);
	const int n = 20;
	for (int i = 0; i < nn; ++i)
	{
		for (int j = 0; j < mm; ++j)
		{
			cout << setw(n) << v[i][j];
		}
		cout << endl;
	}
}

template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> &B) {
	Matrix<T> sum(B.nrows(), B.ncols());
	for (int i = 0; i < B.nrows(); ++i)
	{
		for (int j = 0; j < B.ncols(); ++j)
		{
			sum[i][j] = v[i][j] + B[i][j];
		}
	}
	return sum;
}
template <class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> &B) {
	Matrix<T> diff(B.nrows(), B.ncols());
	for (int i = 0; i < B.nrows(); ++i)
	{
		for (int j = 0; j < B.ncols(); ++j)
		{
			diff[i][j] = v[i][j] - B[i][j];
		}
	}
	return diff;
}

template <class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> &B) {
	Matrix<T> multip(B.nrows(), B.ncols());
	try {
		if (mm != B.nrows()) {
			throw "Error:";
		} else {
			int i, j, k;
			T temp;
			for (i = 0; i < nn; i++) {
				for (j = 0; j < B.nrows(); j++) {
					temp = (0.0, 0.0);
					for (k = 0; k < mm; k++) {
						temp = temp + v[i][k] * B[k][j];
					}
					multip[i][j] = temp;
				}
			}
		}
	}
	catch (const char* msg) {
		cout << msg << " ";
		cout << "matrix product: size doesn't match: ("
		     << nn << "," << mm << "), ("
		     << B.nrows() << "," << B.ncols() << ")." << endl;
	}
	return multip;
}
template <class T>
Matrix<T> Matrix<T>::transpose() {
	Matrix<T> transposed(mm,nn);
	for (int i = 0; i < mm; ++i)
	{
		for (int j = 0; j < nn; ++j)
		{
			transposed[i][j]=v[j][i];
		}
	}
	return transposed;
}
template <class T>
Matrix<T>::~Matrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

template <class T>
Matrix3d<T>::Matrix3d(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
Matrix3d<T>::Matrix3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T * * [n])
{
	int i, j;
	v[0] = new T*[n * m];
	v[0][0] = new T[n * m * k];
	for (j = 1; j < m; j++) v[0][j] = v[0][j - 1] + k;
	for (i = 1; i < n; i++) {
		v[i] = v[i - 1] + m;
		v[i][0] = v[i - 1][0] + m * k;
		for (j = 1; j < m; j++) v[i][j] = v[i][j - 1] + k;
	}
}

template <class T>
inline T** Matrix3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* const * Matrix3d<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int Matrix3d<T>::dim1() const
{
	return nn;
}

template <class T>
inline int Matrix3d<T>::dim2() const
{
	return mm;
}

template <class T>
inline int Matrix3d<T>::dim3() const
{
	return kk;
}

template <class T>
Matrix3d<T>::~Matrix3d()
{
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}
