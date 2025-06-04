#ifndef POSSYMLINSYSTEM_H
#define POSSYMLINSYSTEM_H
#include "matrix.h"
#include "vector.h"
#include "linearsystem.h"

class PosSymLinSystem: public LinearSystem{
    public:
    PosSymLinSystem() = delete; 
    PosSymLinSystem(Matrix* A, Vector* b);
    //solve using conjugate gradient method
    Vector Solve(Vector x0, int max_iters = 100, double min = 0.1);
    //check if symmetric
    bool Check();
};

#endif