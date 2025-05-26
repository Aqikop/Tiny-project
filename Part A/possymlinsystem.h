#ifndef POSSYMLINSYSTEM_H
#define POSSYMLINSYSTEM_H
#include "matrix.h"
#include "vector.h"
#include "linearsystem.h"

class PosSymLinSystem{
    public:
    //solve using conjugate gradient method
    Vector Solve();
    //check if symmetric
    bool Check();
};

#endif