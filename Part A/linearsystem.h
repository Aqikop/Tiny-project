#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H
#include "matrix.h"
#include "vector.h"

class LinearSystem {
private:
    int mSize;
    Matrix* mpA;
    Vector* mpb;
    LinearSystem(const LinearSystem&); 
public:
    LinearSystem() = delete; 
    LinearSystem(Matrix* A, Vector* b);
    Vector Solve();
};

#endif 
