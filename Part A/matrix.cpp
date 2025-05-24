#include <iostream>
#include <cmath>
#include <assert.h>
#include <vector>
using namespace std;

class Matrix{
    private:
        int mNumRows;
        int mNumCols;
        int **mData;
    public:
        //Destructor: Free all the memory
        Matrix(const Matrix&);
        // Constructor accepts two integers
        Matrix(const int, const int);
        //Destructor: Free all the memory
        ~Matrix();
        //# of Rows Getter + # of Columns Getter
        int getNumRows() const;
        int getNumCols() const;
        //Overloadding Operator: ()
        int &operator()(int, int);
        // Unary Operator: +, -
        Matrix operator+() const;
        Matrix operator-() const;
        // Binary Operator: +, -
        Matrix operator+(const Matrix&) const;
        Matrix operator-(const Matrix&) const;
        // Multiplication (with Matrix, Scalar, Vector)
        Matrix operator*(const Matrix&) const;
        Matrix operator*(double) const;
        Matrix Matrix::operator*(const int* ) const;
        //Determinant
        double determinant() const;
        //Inverse
        Matrix inverse() const;
        // Transpose
        Matrix tranpose() const;
        //Pseudo-inverse (Moore-Penrose)
        Matrix pseudo_inverse() const;
};
// Copy Constructor
Matrix::Matrix (const Matrix& matrix){
    mNumRows = matrix.mNumRows;
    mNumCols = matrix.mNumCols;
    mData = new int*[mNumRows];
    for (int i = 0; i < mNumRows; i++){
        mData[i] = new int[mNumCols];
        for (int j = 0; j < mNumCols; j++){
            mData[i][j] = matrix.mData[i][j];
        }
    }
}
 // Constructor accepts two integers
Matrix::Matrix(const int a, const int b){
    mNumRows = a;
    mNumCols = b;
    mData = new int*[mNumRows];
    for (int i = 0; i <mNumRows; i++){
        mData[i] = new int[mNumCols];
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
int &Matrix::operator()(int row, int col){
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
    Matrix mul_matrix(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < matrix.mNumCols; j++){
            for (int k = 0; k < mNumCols; k++){
                mData[i][j] += mData[i][k] * matrix.mData[k][j];
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
Matrix Matrix::operator*(const int* vec) const {
    assert(mNumCols == 1); 
    Matrix vec_mul_matrix(mNumRows, 1);
    for (int i = 0; i < mNumRows; i++) {
        vec_mul_matrix.mData[i][0] = 0;
        for (int j = 0; j < mNumCols; j++) {
            vec_mul_matrix.mData[i][0] += mData[i][j] * vec[j];
        }
    }
    return vec_mul_matrix;
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
            inv.mData[i][j] = static_cast<int>(aug[i][j + n]);
        }
    }
    for (int i = 0; i < n; ++i) delete[] aug[i];
    delete[] aug;
    return inv;
}

//Tranpose
Matrix Matrix::tranpose() const {
    Matrix trans_matrix(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++){
        for (int j = 0; j < mNumCols; j++){
            trans_matrix.mData[j][i] = mData[i][j];
        }
    }
    return trans_matrix;
}
Matrix Matrix::pseudo_inverse() const {
    assert(mNumRows >= mNumCols); 
    Matrix A_trans = this->tranpose(); 
    Matrix AtA = A_trans * (*this);    
    Matrix AtA_inv = AtA.inverse(); 
    Matrix pinv = AtA_inv * A_trans;   
    return pinv;
}


int main(){
    
}