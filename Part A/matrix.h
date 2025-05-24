#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
    int mNumRows;
    int mNumCols;
    int **mData;
public:
    Matrix(const Matrix&);
    Matrix(const int, const int);
    ~Matrix();
    int getNumRows() const;
    int getNumCols() const;
    int &operator()(int, int);

    Matrix operator+() const;
    Matrix operator-() const;
    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix operator*(double) const;
    Matrix operator*(const int*) const;
    double determinant() const;
    Matrix inverse() const;
    Matrix tranpose() const;
    Matrix pseudo_inverse() const;
};

#endif 