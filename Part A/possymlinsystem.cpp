#include <iostream>
#include "matrix.h"
#include "vector.h"
#include "linearsystem.h"
#include "possymlinsystem.h"
#include <cassert>
#include <math.h>
using namespace std;

PosSymLinSystem::PosSymLinSystem(Matrix* A, Vector* b) : LinearSystem(A, b) {
    mSize = A->getNumRows();
    assert(A->getNumRows() == A->getNumCols());
}
Vector PosSymLinSystem::Solve(Vector x0, int max_iters, double min) {
    //initiate some stuff
    Vector r0, p0, x0_temp, r0_temp;
    double r0_dot, alpha, beta;
    x0_temp = x0;
    r0 = (*mpb) - ((*mpA)*x0);
    p0 = r0;
    //sol temp
    int iters = 0;
    while(iters < max_iters){
        //saving old residue
        r0_temp = r0;
        //step size
        r0_dot = r0*r0;  //r0.T * r0
        alpha = r0_dot / (((*mpA)*p0)*p0); // r0_dot / p0.T * mpA * p0
        //update solution
        x0 = x0_temp + alpha*p0; // x0_temp + alpha * p0
        x0_temp = x0;
        //update residue
        r0 = (*mpb) - (*mpA)*x0;
        iters += 1;
        //check stop
        //compute if |r0| < epsilon (very small)
        //r0.length < min => break;
        if(r0.magnitude() < min){break;}
        else{
            beta = (r0*r0) / (r0_temp*r0_temp);
        p0 = r0 + beta*p0;
        }
    }
    return x0;
}

bool PosSymLinSystem::Check() {
    if(mpA->getNumRows() == mpA->getNumCols()){
        double temp = mpA->getNumCols();
        for(int i = 0; i < temp; i++){
            for(int j = i + 1; j < temp; j++){
                if(mpA->getData(i,j) != mpA->getData(j,i)){return false;}
            }
        }
        return true;
    }
    else{return false;}
}
//test
// int main()
// {
//     Matrix m(3,3);
//     m(1,1) = 2;
//     m(1,2) = -1;
//     m(1,3) = 0;
//     m(2,1) = -1;
//     m(2,2) = 2;
//     m(2,3) = -1;
//     m(3,1) = 0;
//     m(3,2) = -1;
//     m(3,3) = 2;
//     Vector k(3);
//     k[0] = 1;
//     k[1] = 0;
//     k[2] = 1;
//     PosSymLinSystem A(&m, &k);
//     Vector p(3);
//     p[0] = 4;
//     p[1] = -1;
//     p[2] = 3;
//     Vector x = A.Solve(p);
//     cout <<"Test 1:";
//     cout << "\n"<<x[0];
//     cout << "\n"<<x[1];
//     cout << "\n"<<x[2];
// }
