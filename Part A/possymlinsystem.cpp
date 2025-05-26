#include <iostream>
#include "matrix.h"
#include "vector.h"
#include "linearsystem.h"
#include "possymlinsystem.h"
#include <cassert>
#include <math.h>
using namespace std;

Vector PosSymLinSystem::Solve(Vector x0 = Vector(mpb->mSize)) {
    //not implement yet
    Vector r0;
    r0 = mpb - mpA*x0
    //
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