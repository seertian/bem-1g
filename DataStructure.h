#ifndef _DataType_H_
#define _DataType_H_

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <math.h>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

using namespace std;
#define EPS 1.0e-7
const double pi = acos(-1);
// macro-like inline functions

template<class T>
inline T SQR(const T a) {return a * a;}

template<class T>
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
{return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
{return b > a ? float(b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
{return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
{return b < a ? float(b) : (a);}

template<class T>
inline T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
{T dum = a; a = b; b = dum;}

// Vector and Matrix Classes

#ifdef _USESTDVECTOR_
#define Vector vector
#else

template <class T>
class Vector {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	Vector();
	explicit Vector(int n);		// Zero-based array
	Vector(int n, const T &a);	//initialize to constant value
	Vector(int n, const T *a);	// Initialize to array
	Vector(const Vector &rhs);	// Copy constructor
	Vector & operator=(const Vector &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	~Vector();
};

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

#endif //ifdef _USESTDVECTOR_

template <class T>
class Matrix {
private:
	int nn;
	int mm;
	T **v;
public:
	Matrix();
	Matrix(int n, int m);			// Zero-based array
	Matrix(int n, int m, const T &a);	//Initialize to constant
	Matrix(int n, int m, const T *a);	// Initialize to array
	Matrix(const Matrix &rhs);		// Copy constructor
	Matrix & operator=(const Matrix &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	//inline Vector* operator[](const int i);
	inline const T* operator[](const int i) const;
	//inline const Vector* Matrix<T>::operator[](const int i) const
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	void print();
	Matrix operator+(Matrix &);
	Matrix operator-(Matrix &);
	Matrix operator*(Matrix &);
	Matrix transpose();
	Matrix normalize_by_row();
	Matrix normalize_by_col();
	~Matrix();
};

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

// template <class T> 
// Vector<T> Matrix<T>::operator[](int num) {
// 	Vector<T> tmp(mm);
// 	for (int  = 0; i < count; ++i)
// 	{
// 		tmp[i]=v[i][num];
// 	}
// 	return tmp;
// }


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
	Matrix<T> transposed(mm, nn);
	for (int i = 0; i < mm; ++i)
	{
		for (int j = 0; j < nn; ++j)
		{
			transposed[i][j] = v[j][i];
		}
	}
	return transposed;
}

template <class T>
Matrix<T> Matrix<T>::normalize_by_row() {
	Matrix<T> normalized(nn, mm);
	T tmp[nn];
	for (int i = 0; i < nn; ++i) {
		tmp[i] = (0.0, 0.0);
		for (int j = 0; j < mm; ++j) {
			tmp[i] = tmp[i] + this->v[i][j] * this->m_matrix[i][j];
		}
	}
	for (int i = 0; i < nn; ++i) {
		for (int j = 0; j < mm; ++j) {
			normalized(i, j) = this->v[i][j] / sqrt(tmp[i]);
		}
	}
	return normalized;
}
template <class T>
Matrix<T> Matrix<T>::normalize_by_col() {
	Matrix normalized(nn, mm, 0.0);
	T tmp[mm];
	for (int i = 0; i < mm; ++i) {
		tmp[i] = (0.0, 0.0);
		for (int j = 0; j < nn; ++j) {
			tmp[i] = tmp[i] + this->v[j][i] * this->m_matrix[j][i];
		}
	}
	for (int i = 0; i < mm; ++i) {
		for (int j = 0; j < nn; ++j) {
			normalized(j, i) = this->v[j][i] / sqrt(tmp[i]);
		}
	}
	return normalized;
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
class Matrix3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	Matrix3d();
	Matrix3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~Matrix3d();
};

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
// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;
typedef double Double;
typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const Vector<Int> VecInt_I;
typedef Vector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const Vector<Uint> VecUint_I;
typedef Vector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const Vector<Llong> VecLlong_I;
typedef Vector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const Vector<Ullong> VecUllong_I;
typedef Vector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const Vector<Char> VecChar_I;
typedef Vector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const Vector<Char*> VecCharp_I;
typedef Vector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const Vector<Uchar> VecUchar_I;
typedef Vector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const Vector<Doub> VecDoub_I;
typedef Vector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const Vector<Doub*> VecDoubp_I;
typedef Vector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const Vector<Complex> VecComplex_I;
typedef Vector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const Vector<Bool> VecBool_I;
typedef Vector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const Matrix<Int> MatInt_I;
typedef Matrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const Matrix<Uint> MatUint_I;
typedef Matrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const Matrix<Llong> MatLlong_I;
typedef Matrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const Matrix<Ullong> MatUllong_I;
typedef Matrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const Matrix<Char> MatChar_I;
typedef Matrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const Matrix<Uchar> MatUchar_I;
typedef Matrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const Matrix<Doub> MatDoub_I;
typedef Matrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const Matrix<Complex> MatComplex_I;
typedef Matrix<Complex> MatComplex, MatComplex_O, MatComplex_IO;

typedef const Matrix<Bool> MatBool_I;
typedef Matrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const Matrix3d<Complex> Mat3DComplex_I;
typedef Matrix3d<Complex> Mat3DComplex, Mat3DComplex_O, Mat3DComplex_IO;
typedef const Matrix3d<Doub> Mat3DDoub_I;
typedef Matrix3d<Complex> Mat3DComplex, Mat3DComplex_O, Mat3DComplex_IO;


#endif /* _DataType_H_ */

