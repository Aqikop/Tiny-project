#include <iostream>
#include "matrix.h"
#include "vector.h"
#include "linearsystem.h"
#include <vector>
#include <cassert>
#include <math.h>
using namespace std;

LinearSystem::LinearSystem(Matrix* A, Vector* b) : mpA(A), mpb(b) {
    mSize = A->getNumRows();
    //assert(A->getNumRows() == A->getNumCols());//may need remove
}
//Solve Method: (Reference: Copilot)
Vector LinearSystem::Solve() {
    if(mpA->getNumRows() == mpA->getNumCols()){
        int n = mSize;
        vector<vector<double>> a(n, vector<double>(n));
        Vector b = *mpb;
        for (int i = 0; i < n; ++i) {
            b[i] = (*mpb)[i];
            for (int j = 0; j < n; ++j) {
                a[i][j] = (*mpA)(i+1, j+1); // 1-based indexing in Matrix
            }
        }
        // Gaussian elimination with partial pivoting
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i+1; j < n; ++j) {
                if (fabs(a[j][i]) > fabs(a[pivot][i])) pivot = j;
            }
            assert(fabs(a[pivot][i]) > 1e-12);
            if (pivot != i) {
                swap(a[i], a[pivot]);
                swap(b[i], b[pivot]);
            }
            for (int j = i+1; j < n; ++j) {
                double factor = a[j][i] / a[i][i];
                for (int k = i; k < n; ++k) {
                    a[j][k] -= factor * a[i][k];
                }
                b[j] -= factor * b[i];
            }
        }
        // Back substitution
        Vector res(n);
        for (int i = n-1; i >= 0; --i) {
            res[i] = b[i];
            for (int j = i+1; j < n; ++j) {
                res[i] -= a[i][j] * res[j];
            }
            res[i] /= a[i][i];
        }
        return res;
    }
    else if(mpA->getNumRows() > mpA->getNumCols()){
        //for overdetermined
        return ((*mpA).pseudo_inverse() * (*mpb));
    }
    else{
        //for underdetermined
        return ((*mpA).pseudo_inverse() * (*mpb));
    }
}

Vector LinearSystem::Solve(double lambda) {
    if(mpA->getNumRows() > mpA->getNumCols()){
        //for overdetermined
        return ((*mpA).pseudo_inverse(lambda) * (*mpb));
    }
    else{
        //for underdetermined
        return ((*mpA).pseudo_inverse(lambda) * (*mpb));
    }
}
//test
int main()
{
    cout <<"test";
    Matrix m(2,3);
    m(1,1) = 1;
    m(1,2) = 2;
    m(1,3) = -1;
    m(2,1) = 3;
    m(2,2) = -1;
    m(2,3) = 2;
    Vector k(2);
    k[0] = 1;
    k[1] = 5;
    LinearSystem A(&m, &k);
    Vector x = A.Solve(1);
    cout << "\n"<<x[0];
    cout << "\n"<<x[1];
    cout << "\n"<<x[2];
}