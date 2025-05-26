#include <iostream>
#include "matrix.h"
#include "vector.h"
#include <vector>
#include <cassert>
#include <math.h>
using namespace std;

class LinearSystem {
    protected:
        int mSize;
        Matrix* mpA;
        Vector* mpb;
        // Prevent copy constructor
        LinearSystem(const LinearSystem&);
    public:
        // No default constructor
        LinearSystem() = delete;
        // Constructor requires matrix and vector
        LinearSystem(Matrix* A, Vector* b);
        // Solve method 
        virtual Vector Solve(); 
};

LinearSystem::LinearSystem(Matrix* A, Vector* b) : mpA(A), mpb(b) {
    mSize = A->getNumRows();
    assert(A->getNumRows() == A->getNumCols());
}
//Solve Method: (Reference: Copilot)
Vector LinearSystem::Solve() { 
    int n = mSize;
    vector<std::vector<double>> a(n, vector<double>(n));
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
