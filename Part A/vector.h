#ifndef VECTOR_H
#define VECTOR_H

class Vector{
    private:
        int mSize;
        double* mData;
    public:
        // Constructor and Destructor
        Vector();
        Vector(int size);
        Vector(const Vector& other);
        ~Vector();

        Vector& operator=(const Vector& other); // Assignment operator

        // Unary operators
        Vector operator-() const;
        Vector operator+() const;
        Vector operator++();
        Vector operator--();

        // Binary operators
        Vector operator+(const Vector& other) const;
        Vector operator-(const Vector& other) const;
        Vector operator*(double scalar) const;
        double operator*(const Vector& other) const; // Dot product(v₁·v₂ = a×d + b×e + c×f )

        double& operator[](int index);
        double& operator()(int index);
        double operator[](int index) const {
            return mData[index];
        }

        int get_size() const{
            return mSize;
        };
};
Vector operator*(double scalar, const Vector& other); // Scalar multiplication (v₁×a = a×v₁ = [a×b, a×c, a×d])

#endif
