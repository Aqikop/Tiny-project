#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include "vector.h"

class Matrix {
private:
    int mNumRows;
    int mNumCols;
    double **mData;
public:
    Matrix(const Matrix&);
    Matrix(const int, const int);
    ~Matrix();
    int getNumRows() const;
    int getNumCols() const;
    int get(int row, int col) const { return mData[row][col];}
    double &operator()(int, int);

    Matrix operator+() const;
    Matrix operator-() const;
    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix operator*(double) const;
    //Matrix operator*(const int*) const;
    double determinant() const;
    Matrix inverse() const;
    Matrix tranpose() const;
    Matrix pseudo_inverse() const;
};
Vector operator*(const Matrix& m, const Vector& v); //row mul
Vector operator*(const Vector& v, const Matrix& m); //col mul

#endif 