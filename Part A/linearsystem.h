#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H
#include "matrix.h"
#include "vector.h"

class LinearSystem {
protected:
    int mSize;
    Matrix* mpA;
    Vector* mpb;
    LinearSystem(const LinearSystem&); 
public:
    LinearSystem() = delete; 
    LinearSystem(Matrix* A, Vector* b);
    virtual Vector Solve();
};

#endif 
