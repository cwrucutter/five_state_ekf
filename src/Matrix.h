#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#endif

#include <iostream>
#include <vector>

template<class T>
class Matrix {
private:
	unsigned int rows_, cols_;
	T covariance(const std::vector<T>&, const T&, const std::vector<T>&, const T&);
public:
	Matrix();
	Matrix(int, int);
	Matrix(int, int, const T&);
	Matrix(const Matrix<T>&);
	virtual ~Matrix() {}

	Matrix<T> diagonalMatrix(int, int, T&);
	Matrix<T> diagonalMatrix(int, int, T);
	Matrix<T> identity(int);

	Matrix<T>& operator= (const Matrix<T>&);
	Matrix<T> operator+ (const Matrix<T>&);
	Matrix<T>& operator+= (const Matrix<T>&);
	Matrix<T> operator- (const Matrix<T>&);
	Matrix<T>& operator-= (const Matrix<T>&);
	Matrix<T> operator* (const Matrix<T>&);
	Matrix<T>& operator*= (const Matrix<T>&);
	bool operator==(const Matrix<T>&);
	bool equals(int, int, const std::vector<std::vector<T> >&);
	Matrix<T> transpose(void);
	Matrix<T> covariance(void);
	Matrix<T> inverse(void);

	Matrix<T> operator+ (const T&);
	Matrix<T> operator- (const T&);
	Matrix<T> operator* (const T&);
	Matrix<T> operator/ (const T&);

	std::vector<T> operator* (const std::vector<T>&);

	T& operator() (const unsigned int&, const unsigned int&);
	const T& operator() (const unsigned int&, const unsigned int&) const;

	T& operator() (const int&, const int&);
	const T& operator() (const int&, const int&) const;

	unsigned int getRows(void) const { return this->rows_; }
	unsigned int getColumns(void) const { return this->cols_; }

	void print(void);
	void printToFile(std::ofstream&);
protected:
	std::vector<std::vector<T> > matrix_;

};

template<class T>
Matrix<T>::Matrix() {
	this->rows_ = this->cols_ = 2;
	this->matrix_.resize(2);
	for (int i = 0; i < 2; i++) { this->matrix_[i].resize(2); }
}

template<class T>
Matrix<T>::Matrix(int rows, int cols) {
	this->matrix_.resize(rows);
	for (int i = 0; i < rows; i++) { this->matrix_[i].resize(cols); }
	this->rows_ = rows;
	this->cols_ = cols;
}

template<class T>
Matrix<T>::Matrix(int rows, int cols, const T& initVal) {
	this->matrix_.resize(rows);
	for (int i = 0; i < rows; i++) { this->matrix_[i].resize(cols, initVal); }
	this->rows_ = rows;
	this->cols_ = cols;
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& m) {
	this->matrix_ = m.matrix_;
	this->rows_ = m.getRows();
	this->cols_ = m.getColumns();
}

template<class T>
Matrix<T> Matrix<T>::identity(int size) {
	Matrix<T> result( size, size, 0.0 );
	for (unsigned int i = 0; i < result.getRows(); i++) {
		result(i, i) = 1.0;
	}
	return result;
}

template<class T>
Matrix<T> Matrix<T>::diagonalMatrix(int rows, int cols, T& diagVal) {
	Matrix<T> result( rows, cols, 0.0 );
	for (unsigned int i = 0; i < result.getRows(); i++) { result(i, i) = diagVal; }
	return result;
}

template<class T>
Matrix<T> Matrix<T>::diagonalMatrix(int rows, int cols, T diagVal) {
	Matrix<T> result( rows, cols, 0.0 );
	for (unsigned int i = 0; i < result.getRows(); i++) { result(i, i) = diagVal; }
	return result;
}

template<class T>
T& Matrix<T>::operator() (const unsigned int& i, const unsigned int& j) { return this->matrix_[i][j]; }

template<class T>
const T& Matrix<T>::operator() (const unsigned int& i, const unsigned int& j) const { return this->matrix_[i][j]; }

template<class T>
T& Matrix<T>::operator() (const int& i, const int& j) { return this->matrix_[i][j]; }

template<class T>
const T& Matrix<T>::operator() (const int& i, const int& j) const { return this->matrix_[i][j]; }

template<class T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& rhs) {
	if (this == &rhs) { return *this; } //why bother adding stuff to myself?
	unsigned int rows = rhs.getRows();
	unsigned int cols = rhs.getColumns();

	this->matrix_.resize(rows);
	for (unsigned int i = 0; i < rows; i++) { this->matrix_[i].resize(cols); }

	for (unsigned int i = 0; i < rows; i++) {
		for (unsigned int j = 0; j < cols; j++) {
			this->matrix_[i][j] = rhs(i, j);
		}
	}
	this->rows_ = rows;
	this->cols_ = cols;
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& rhs) {
	Matrix<T> result( (int)this->getRows(), (int)this->getColumns(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result(i, j) = this->matrix_[i][j] + rhs(i, j);
		}
	}
	return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator+= (const Matrix<T>& rhs) {
	for (unsigned int i = 0; i < rhs.getRows(); i++) {
		for (unsigned int j = 0; j < rhs.getColumns(); j++) {
			this->matrix_[i][j] += rhs(i, j);
		}
	}
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator- (const Matrix<T>& rhs) {
	Matrix<T> result( (int)this->getRows(), (int)this->getColumns(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result(i, j) = this->matrix_[i][j] - rhs(i, j);
		}
	}
	return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator-= (const Matrix<T>& rhs) {
	for (unsigned int i = 0; i < rhs.getRows(); i++) {
		for (unsigned int j = 0; j <= rhs.getColumns(); j++) {
			this->matrix_[i][j] -= rhs(i, j);
		}
	}
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& rhs) {
	Matrix<T> result( (int)this->getRows(), (int)rhs.getColumns(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < rhs.getColumns(); j++) {
			for (unsigned int k = 0; k < this->getColumns(); k++) {
				result(i, j) += this->matrix_[i][k] * rhs(k, j);
			}
		}
	}
	return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator*= (const Matrix<T>& rhs) {
	Matrix<T> result = (*this) * rhs;
	*this = result;
	return *this;
}

template<class T>
bool Matrix<T>::operator== (const Matrix<T>& rhs) {
	if (this->getRows() != rhs.getRows() || this->getColumns() != rhs.getColumns()) { return false; }
	for (unsigned int i = 0; i < rhs.getRows(); i++) {
		for (unsigned int j = 0; j < rhs.getColumns(); j++) { if ((*this)(i, j) != rhs(i, j)) { return false; } }
	}
	return true;
}

template<class T>
bool Matrix<T>::equals(int rows, int cols, const std::vector<std::vector<T> >& vals) {
	if (rows != this->getRows() || cols != this->getColumns()) { return false; }
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if ((*this)(i, j) != vals[i][j]) { return false; }
		}
	}
	return true;
}

template<class T>
Matrix<T> Matrix<T>::transpose(void) {
	Matrix<T> result( (int)this->getColumns(), (int)this->getRows(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result(j, i) = this->matrix_[i][j];
		}
	}
	return result;
}

template<class T>
T Matrix<T>::covariance(const std::vector<T>& x, const T& xMean, const std::vector<T>& y, const T& yMean) {
	T covar = 0.0;
	for (unsigned int i = 0; i < x.size(); i++) {
		covar += (x[i] - xMean)*(y[i] - yMean);
	}
	return covar;
}

template<class T>
Matrix<T> Matrix<T>::covariance(void) {
	T mean = 0.0, covar = 0.0;
	std::vector<T> means( this->getRows(), 0.0 );
	for (unsigned int i = 0; i <= this.getRows(); i++) {
		for (unsigned int j = 0; j < this.getColumns(); j++) {
			mean += (*this)(i, j);
		}
		means.push_back(mean / this->getColumns());
	}

	Matrix<T> result( this->getRows(), this->getRows(), 0.0 );
	for (unsigned int i = 0; i < this.getRows(); i++) {
		for (unsigned int j = i; j < this.getRows(); j++) {
			covar = this->covariance(this->matrix_[i], means[i], this->matrix_[j], means[j]);
			this(i, j) = covar;
			this(j, i) = covar;
		}
	}
	return mean;
}

template<class T>
Matrix<T> Matrix<T>::inverse(void) {
	//use Gauss-Jordan Elimination
	Matrix<T> result;
	if (this->getRows() != this->getColumns()) { return *this; }
	if (this->getRows() == this->getColumns() && this->getRows() == 1) {
		return Matrix < T >(1, 1, 1.0 / (*this)(0, 0));
	}
	else if (this->getRows() == this->getColumns() && this->getRows() == 2) {
		T determinant = (*this)(0, 0)*(*this)(1, 1) - (*this)(0, 1)*(*this)(1, 0);
		result = Matrix < T >((int)this->getRows(), (int)this->getColumns());
		result(0, 0) = (*this)(1, 1) / determinant;
		result(1, 1) = (*this)(0, 0) / determinant;
		result(0, 1) = -1 * (*this)(0, 1) / determinant;
		result(1, 0) = -1 * (*this)(1, 0) / determinant;
		return result;
	}
	else {
		return *this;
	}
}

template<class T>
Matrix<T> Matrix<T>::operator+ (const T& rhs) {
	Matrix<T> result( this->getRows(), this->getColumns(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result(i, j) = this->matrix_[i][j] + rhs;
		}
	}
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator- (const T& rhs) {
	Matrix<T> result( this->getRows(), this->getColumns(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result(i, j) = this->matrix_[i][j] - rhs;
		}
	}
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator* (const T& rhs) {
	Matrix<T> result( this->getRows(), this->getColumns(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result(i, j) = this->matrix_[i][j] * rhs;
		}
	}
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator/ (const T& rhs) {
	Matrix<T> result( this->getRows(), this->getColumns(), 0.0 );
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result(i, j) = this->matrix_[i][j] / rhs;
		}
	}
	return result;
}

template<class T>
std::vector<T> Matrix<T>::operator* (const std::vector<T>& rhs) {
	std::vector<T> result;
	result.resize(this->getRows(), 0.0);
	for (unsigned int i = 0; i < this->getRows(); i++) {
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			result[i] += this->matrix_[i][j] * rhs[j];
		}
	}
	return result;
}

template<class T>
void Matrix<T>::print(void) {
	for (unsigned int i = 0; i < this->getRows(); i++) {
		std::cout << "[";
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			std::cout << (*this)(i, j);
			if (j < this->getColumns() - 1) { std::cout << ", "; }
		}
		std::cout << "]" << std::endl;
	}
}

template<class T>
void Matrix<T>::printToFile(std::ofstream& file) {
	for (unsigned int i = 0; i < this->getRows(); i++) {
		file << "[";
		for (unsigned int j = 0; j < this->getColumns(); j++) {
			file << (*this)(i, j);
			if (j < this->getColumns() - 1) { file << ", "; }
		}
		file << "]" << std::endl;
	}
}