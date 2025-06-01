#include <iostream>
#include <cmath>
#include <assert.h>
#include <vector>
#include "matrix.h"
#include "vector.h"
using namespace std;

// Copy Constructor
Matrix::Matrix (const Matrix& matrix){
    mNumRows = matrix.mNumRows;
    mNumCols = matrix.mNumCols;
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; i++){
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; j++){
            mData[i][j] = matrix.mData[i][j];
        }
    }
}
 // Constructor accepts two integers
Matrix::Matrix(const int a, const int b){
    mNumRows = a;
    mNumCols = b;
    mData = new double*[mNumRows];
    for (int i = 0; i <mNumRows; i++){
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; j++){
            mData[i][j] = 0;
        }
    }
}
//Destructor: Free all the memory
Matrix::~Matrix(){
    for (int i = 0; i < mNumRows; i++){
        delete[] mData[i]; 
    }
    delete[] mData; 
}
//# of Rows Getter + # of Columns Getter
int Matrix::getNumRows() const{
    return mNumRows;
}
int Matrix::getNumCols() const{
    return mNumCols;
}
// Overloadding Operator: ()
double &Matrix::operator()(int row, int col){
    return mData[row - 1][col - 1];
}
// Overloadding Operator: =


//Unary Operator:
Matrix Matrix::operator+() const{
    return *this;
}
Matrix Matrix::operator-() const {
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < mNumCols; j++){
            mData[i][j] = - mData[i][j];
        }
    }
    return *this;
}
// Addition +
Matrix Matrix::operator+(const Matrix& matrix) const{
    assert(mNumRows == matrix.mNumRows && mNumCols == matrix.mNumCols);
    Matrix sum_matrix(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < mNumCols; j++){
            sum_matrix.mData[i][j] = mData[i][j] + matrix.mData[i][j];
        }
    }
    return sum_matrix;
}
// Subtraction -
Matrix Matrix::operator-(const Matrix & matrix) const{
    assert(mNumRows == matrix.mNumRows && mNumCols == matrix.mNumCols);
    Matrix diff_matrix(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < mNumCols; j++){
            diff_matrix.mData[i][j] = mData[i][j] - matrix.mData[i][j];
        }
    }
    return diff_matrix;
}

//Multiplication *: With matrix, scalar and vector;
Matrix Matrix::operator*(const Matrix & matrix) const{
    assert(mNumCols == matrix.mNumRows);
    Matrix mul_matrix(mNumRows, matrix.mNumCols);
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < matrix.mNumCols; j++){
            mul_matrix.mData[i][j] = 0;
            for (int k = 0; k < mNumCols; k++){
                mul_matrix.mData[i][j] += mData[i][k] * matrix.mData[k][j];
            }
        }
    }
    return mul_matrix;
}
Matrix Matrix::operator*(double scalar) const {
    Matrix scalar_mul_matrix(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < mNumCols; j++){
            scalar_mul_matrix.mData[i][j] = mData[i][j] * scalar;
        }
    }
    return scalar_mul_matrix;
}

// Row multiplication( matrix * vector )
Vector operator*(const Matrix& m, const Vector& v){
    assert(m.getNumCols() == v.get_size());
    Vector result(m.getNumRows());
    for(int i = 0; i < m.getNumRows(); i++){
        double sum = 0;
        for(int j = 0; j < m.getNumCols(); j++){
            sum += m.getData(i,j) * v[j];
        }
        result[i] = sum;
    }
    return result;
}

// Col multiplication (vector * matrix)
Vector operator*(const Vector& v, const Matrix& m){
    assert(v.get_size() == m.getNumRows());
    Vector result(m.getNumCols());
    for(int i = 0; i < m.getNumCols(); i++){
        double sum = 0;
        for(int j = 0; j < m.getNumRows(); j++){
            sum += v[j] * m.getData(j,i);
        }
        result[i] = sum;
    }
    return result;
}

//Determinant (Reference: GPT)
double Matrix::determinant() const {
    if (mNumRows != mNumCols) return 0.0;
    int n = mNumRows;
    double det = 1.0;
    // Create a copy of the matrix as double
    double **temp = new double*[n];
    for (int i = 0; i < n; ++i) {
        temp[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            temp[i][j] = mData[i][j];
        }
    }
    for (int i = 0; i < n; ++i) {
        // Find pivot
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(temp[j][i]) > fabs(temp[pivot][i])) pivot = j;
        }
        if (fabs(temp[pivot][i]) < 1e-12) {
            // Singular matrix
            for (int k = 0; k < n; ++k) delete[] temp[k];
            delete[] temp;
            return 0.0;
        }
        if (pivot != i) {
            std::swap(temp[i], temp[pivot]);
            det = -det;
        }
        det *= temp[i][i];
        for (int j = i + 1; j < n; ++j) {
            temp[i][j] /= temp[i][i];
        }
        for (int j = i + 1; j < n; ++j) {
            for (int k = i + 1; k < n; ++k) {
                temp[j][k] -= temp[j][i] * temp[i][k];
            }
        }
    }
    for (int k = 0; k < n; ++k) delete[] temp[k];
    delete[] temp;
    return det;
}

// Inverse (Reference: GPT)
Matrix Matrix::inverse() const {
    assert(mNumRows == mNumCols);
    int n = mNumRows;
    Matrix inv(n, n);
    double **aug = new double*[n];
    for (int i = 0; i < n; ++i) {
        aug[i] = new double[2 * n];
        for (int j = 0; j < n; ++j) {
            aug[i][j] = mData[i][j];
            aug[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(aug[j][i]) > fabs(aug[pivot][i])) pivot = j;
        }
        assert(fabs(aug[pivot][i]) > 1e-12);
        if (pivot != i) {
            std::swap(aug[i], aug[pivot]);
        }
        double div = aug[i][i];
        for (int j = 0; j < 2 * n; ++j) aug[i][j] /= div;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double factor = aug[j][i];
                for (int k = 0; k < 2 * n; ++k) {
                    aug[j][k] -= factor * aug[i][k];
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inv.mData[i][j] = (aug[i][j + n]);
        }
    }
    for (int i = 0; i < n; ++i) delete[] aug[i];
    delete[] aug;
    return inv;
}

//Tranpose
Matrix Matrix::tranpose() const {
    Matrix trans_matrix(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < mNumCols; j++){
            trans_matrix.mData[j][i] = mData[i][j];
        }
    }
    return trans_matrix;
}
// Pseudo inverse (Ref: chat GPT)
Matrix Matrix::pseudo_inverse(double lamba = 0) const {
    if (mNumRows >= mNumCols) {
        // Tall or square matrix: A⁺ = (AᵀA)⁻¹Aᵀ
        Matrix A_trans = this->tranpose();
        Matrix AtA = A_trans * (*this);
        //creat a identity matrix with same size of mNumcols x mNumcols
        Matrix I(mNumCols, mNumCols);
        for(int i = 1; i <= mNumCols; i++){
            I(i,i) = 1;
        }
        //calculate the psuedo inverse with Tikhonov regularization
        Matrix AtA = AtA + (I)*(lambda);
        
        if (fabs(AtA.determinant()) < 1e-12) {
            std::cerr << "Matrix is rank deficient, pseudo-inverse cannot be computed\n";
            return Matrix(mNumCols, mNumRows);
        }
        Matrix AtA_inv = AtA.inverse();
        return AtA_inv * A_trans;
    } 
    else {
        // Wide matrix: A⁺ = Aᵀ(AAᵀ)⁻¹
        Matrix A_trans = this->tranpose();
        Matrix AAt = (*this) * A_trans;
        //creat a identity matrix with same size of mNumcols x mNumcols
        Matrix I(mNumCols, mNumCols);
                for(int i = 1; i <= mNumCols; i++){
            I(i,i) = 1;
        }
        //calculate the psuedo inverse with Tikhonov regularization
        Matrix AtA = AtA + (I)*(lambda);
        
        if (fabs(AAt.determinant()) < 1e-12) {
            std::cerr << "Matrix is rank deficient, pseudo-inverse cannot be computed\n";
            return Matrix(mNumCols, mNumRows);
        }
        Matrix AAt_inv = AAt.inverse();
        return A_trans * AAt_inv;
    }
}


//Testing code provided by chatgpt, uncomment to test the each function
int main(){
    // Test constructor and getNumRows/getNumCols
    Matrix m1(2, 3);
    std::cout << "Rows: " << m1.getNumRows() << ", Cols: " << m1.getNumCols() << std::endl;

    // Test operator()
    m1(1, 1) = 5;
    m1(1, 2) = 6;
    m1(2, 1) = 7;
    std::cout << "m1(1,1): " << m1(1, 1) << ", m1(2,1): " << m1(2, 1) << std::endl;

    // Test copy constructor
    Matrix m2 = m1;
    std::cout << "Copy m2(1,1): " << m2(1, 1) << std::endl;

    // Test unary +
    Matrix m3 = +m1;
    std::cout << "Unary + m3(1,1): " << m3(1, 1) << std::endl;

    // Test unary -
    Matrix m4 = -m1;
    std::cout << "Unary - m4(1,1): " << m4(1, 1) << std::endl;

    // Test operator+ (addition)
    Matrix m5(2, 3);
    m5(1, 1) = 1;
    m5(2, 1) = 2;
    Matrix m6 = m1 + m5;
    std::cout << "Addition m6(1,1): " << m6(1, 1) << std::endl;

    // Test operator- (subtraction)
    Matrix m7 = m1 - m5;
    std::cout << "Subtraction m7(1,1): " << m7(1, 1) << std::endl;

    // Test operator* (scalar)
    Matrix m8 = m1 * 2.0;
    std::cout << "Scalar multiply m8(1,1): " << m8(1, 1) << std::endl;

    // Test operator* (matrix multiplication)
    Matrix a(2, 2);
    Matrix b(2, 2);
    a(1, 1) = 1; a(1, 2) = 2; a(2, 1) = 3; a(2, 2) = 4;
    b(1, 1) = 2; b(1, 2) = 0; b(2, 1) = 1; b(2, 2) = 2;
    Matrix c = a * b;
    std::cout << "Matrix multiply c(1,1): " << c(1, 1) << std::endl;

    // Test Matrix * Vector (column vector)
    Matrix mat1(2, 3);
    mat1(1, 1) = 1; mat1(1, 2) = 2; mat1(1, 3) = 3;
    mat1(2, 1) = 4; mat1(2, 2) = 5; mat1(2, 3) = 6;
    Vector vec1(3);
    vec1[0] = 1; vec1[1] = 0; vec1[2] = -1;
    Vector result1 = mat1 * vec1;
    std::cout << "Matrix * Vector result: [" << result1[0] << ", " << result1[1] << "]" << std::endl;
    // Expected: [1*1 + 2*0 + 3*(-1) = -2, 4*1 + 5*0 + 6*(-1) = -2]

    // Test Vector * Matrix (row vector)
    Vector vec2(2);
    vec2[0] = 2; vec2[1] = 3;
    Matrix mat2(2, 3);
    mat2(1, 1) = 1; mat2(1, 2) = 2; mat2(1, 3) = 3;
    mat2(2, 1) = 4; mat2(2, 2) = 5; mat2(2, 3) = 6;
    Vector result2 = vec2 * mat2;
    std::cout << "Vector * Matrix result: [" << result2[0] << ", " << result2[1] << ", " << result2[2] << "]" << std::endl;

    // Test determinant
    Matrix detMat(2, 2);
    detMat(1, 1) = 1; detMat(1, 2) = 2; detMat(2, 1) = 3; detMat(2, 2) = 4;
    std::cout << "Determinant: " << detMat.determinant() << std::endl;

    // Test inverse (check)
    Matrix invMat = detMat.inverse();
    std::cout << "Inverse invMat(1,1): " << invMat(1, 1) << std::endl;

    // Test transpose
    Matrix tMat = detMat.tranpose();
    std::cout << "Transpose tMat(1,2): " << tMat(1, 2) << std::endl;

    // Test pseudo-inverse (check)
    Matrix pinv = detMat.pseudo_inverse();
    std::cout << "Pseudo-inverse pinv(1,1): " << pinv(1, 1) << std::endl;

    return 0;
}