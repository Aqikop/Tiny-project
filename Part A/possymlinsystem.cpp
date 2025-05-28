#include <iostream>
#include "matrix.h"
#include "vector.h"
#include "linearsystem.h"
#include "possymlinsystem.h"
#include <cassert>
#include <math.h>
using namespace std;

Vector PosSymLinSystem::Solve(Vector x0 = Vector(mpb->mSize), int max_iters = 100, double min = 0.1) {
    //not implement yet
    Vector r0, p0, x0_temp, r0_temp;
    double r0_dot, alpha, beta;
    x0_temp = x0;
    r0 = mpb - mpA*x0;
    p0 = r0;
    //sol temp
    int iters = 0;
    while(iters < max_iters){
        //saving old residue
        r0_temp = r0;
        //step size
        r0_dot = r0*r0;  //r0.T * r0
        alpha = r0_dot / ((mpA*p0)*p0); // r0_dot / r0.T * mpA * r0
        //update solution
        x0 = x0_temp + alpha*p0; // x0_temp + alpha * p0
        x0_temp = x0;
        //update residue
        r0 = mpb - mpA*x0;
        iters += 1
        //check stop
        //compute if |r0| < epsilon (very small)
        //r0.length < min => break;
        if(r0.length() < min){break;}
        else{
            beta = (r0*r0) / (r0_temp*r0_temp);
        p0 = r0 + beta*p0;
        }
    }
    return x0;
}

bool PosSymLinSystem::Check() {
    if(mpA->mNumRows == mpA->mNumCols){
        double temp = mpA->mNumRows;
        for(int i = 0; i < mNumRows; i++){
            for(int j = i + 1; j < mNumRows; j++){
                if(mData[i][j] != mData[j][i]){return false}
            }
        }
        return true
    }
    else{return false}
}