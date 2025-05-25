#include <iostream>
#include "vector.h"
using namespace std;

Vector:: Vector() : mSize(0), mData(nullptr) {
}

Vector:: Vector(int size): mSize(size){
    mData = new double[size];
    for(int i = 0; i < size; i++){
        mData[i] = 0.0;
    }
}

Vector:: Vector(const Vector& other): mSize(other.mSize){
    mData = new double[mSize];
    for(int i = 0; i < mSize; i++){
        mData[i] = other.mData[i];
    }
}

Vector:: ~Vector(){
    delete[] mData;
}

// Assignment operator
Vector& Vector :: operator=(const Vector& other){
    if(this == &other){
        return *this;
    }
    mSize = other.mSize;
    mData = new double[mSize];
    for(int i = 0; i < mSize; i++){
        mData[i] = other.mData[i];
    }
    return *this;
}

// Unary operators
Vector Vector:: operator-() const{
    Vector temp(mSize);
    for(int i = 0; i < mSize; i++){
        temp.mData[i] = -mData[i];
    }
    return temp;
}

Vector Vector:: operator+() const{
    Vector temp(mSize);
    for(int i = 0; i < mSize; i++){
        temp.mData[i] = +mData[i];
    }
    return temp;
}

Vector Vector:: operator++(){
    for(int i = 0; i < mSize; i++){
        mData[i]++;
    }
    return *this;
}

Vector Vector:: operator--(){
    for(int i = 0; i < mSize; i++){
        mData[i]--;
    }
    return *this;
}

//Binary operations
Vector Vector:: operator+(const Vector& other) const{
    Vector result(mSize);
    if(mSize != other.mSize){
        cout << "2 vectors have different size, cannot compute";
    }
    for(int i = 0; i < mSize; i++){
        result.mData[i] = mData[i] + other.mData[i];
    }
    return result;
}

Vector Vector:: operator-(const Vector& other) const{
    Vector result(mSize);
    if(mSize != other.mSize){
        cout << "2 vectors have different size, cannot compute";
    }
    for(int i = 0; i < mSize; i++){
        result.mData[i] = mData[i] - other.mData[i];
    }
    return result;
}

Vector Vector:: operator*(double scalar) const{
    Vector result(mSize);
    for(int i = 0; i < mSize; i++){
        result.mData[i] = scalar * mData[i];
    }
    return result;
}

double Vector:: operator*(const Vector& other) const{
    double sum;
    if(mSize != other.mSize){
        cout << "2 vectors have different size, cannot compute";
    }
    for(int i = 0; i < mSize; i++){
        sum += mData[i] * other.mData[i];
    }
    return sum;
}

double& Vector:: operator[](int index){
    if(index < 0 || index > mSize){
        cout << "Index is out of bound";
    }
    return mData[index];
}

double& Vector:: operator()(int index){
    if(index < 0 || index > mSize){
        cout << "Index is out of bound";
    }
    return mData[index - 1]; // vị trí đầu tiên là 0 nhưng mà bth ko biết thì coi vị trí đầu tiên là 1
}

double Vector:: operator[](int index) const{
    return mData[index];
}

Vector operator*(double scalar, const Vector& v2){
    return v2 * scalar;
};


void printVector(const Vector& v, string name) {
    cout << name << " (size " << v.get_size() << "): [";
    for (int i = 0; i < v.get_size(); i++) {
        // Need to cast away constness since operator[] isn't const
        cout << const_cast<Vector&>(v)[i];
        if (i < v.get_size() - 1) cout << ", ";
    }
    cout << "]" << endl;
}

// // Testing code provided by chatgpt
// int main(){
//     // Test default constructor
//     cout << "Testing default constructor:" << endl;
//     Vector v1;
//     cout << "Default vector size: " << v1.get_size() << endl;

//     // Test parameterized constructor
//     cout << "\nTesting parameterized constructor:" << endl;
//     Vector v2(5);
//     printVector(v2, "v2");

//     // Initialize vector with values
//     for (int i = 0; i < v2.get_size(); i++) {
//         v2[i] = i + 1.0;
//     }
//     printVector(v2, "v2 (after setting values)");

//     // Test copy constructor
//     cout << "\nTesting copy constructor:" << endl;
//     Vector v3(v2);
//     printVector(v3, "v3 (copy of v2)");

//     // Modify v2 to ensure deep copy was made
//     v2[0] = 99.0;
//     printVector(v2, "v2 (after modification)");
//     printVector(v3, "v3 (should be unchanged)");


//     cout << "\nTesting assignment operator:" << endl;
//     Vector v4(3);
//     for (int i = 0; i < v4.get_size(); i++) {
//         v4[i] = (i + 1) * 10.0;
//     }
//     printVector(v4, "v4 (original)");

//     v4 = v2;
//     printVector(v4, "v4 (after v4 = v2)");

//     // Self-assignment
//     v4 = v4;
//     printVector(v4, "v4 (after self-assignment)");

//     cout << "\nTesting unary operators:" << endl;
//     Vector v5(3);
//     for (int i = 0; i < v5.get_size(); i++) {
//         v5[i] = i - 1.0;  // -1, 0, 1
//     }
//     printVector(v5, "v5");

//     // Test unary minus
//     Vector v6 = -v5;
//     printVector(v6, "-v5");

//     // Test unary plus
//     Vector v7 = +v5;
//     printVector(v7, "+v5");

//     // Test pre-increment
//     Vector v8 = v5;
//     printVector(++v8, "++v5");

//     // Test pre-decrement
//     Vector v9 = v5;
//     printVector(--v9, "--v5");


//     cout << "\nTesting binary operators:" << endl;
//     Vector a(3);
//     Vector b(3);

//     // Initialize vectors
//     for (int i = 0; i < 3; i++) {
//         a[i] = i + 1.0;       // 1, 2, 3
//         b[i] = (i + 1) * 2.0; // 2, 4, 6
//     }
//     printVector(a, "a");
//     printVector(b, "b");

//     // Test addition
//     Vector c = a + b;
//     printVector(c, "a + b");

//     // Test subtraction
//     Vector d = a - b;
//     printVector(d, "a - b");

//     // Test scalar multiplication (both ways)
//     Vector e = a * 2.0;
//     printVector(e, "a * 2.0");

//     Vector f = 3.0 * b;
//     printVector(f, "3.0 * b");

//     // Test dot product
//     double dot = a * b;
//     cout << "a * b (dot product): " << dot << endl;

//     cout << "\nTesting indexing operators:" << endl;
//     Vector v10(5);
//     for (int i = 0; i < v10.get_size(); i++) {
//         v10[i] = (i + 1) * 1.5;
//     }
//     printVector(v10, "v10");

//     // Test [] operator
//     cout << "v10[0] = " << v10[0] << ", v10[2] = " << v10[2] << endl;

//     // Test () operator (1-based indexing)
//     cout << "v10(1) = " << v10(1) << ", v10(3) = " << v10(3) << endl;

//     // Modify using both operators
//     v10[1] = 100.0;
//     v10(5) = 200.0;
//     printVector(v10, "v10 (after modification)");

//     cout << "\nTesting error handling:" << endl;
//     // Test vectors of different sizes
//     Vector small(2);
//     Vector large(4);
//     small[0] = small[1] = 1.0;
//     for (int i = 0; i < large.get_size(); i++) large[i] = i;

//     cout << "Adding vectors of different sizes:" << endl;
//     Vector result = small + large;  // Should print error message

//     // Test bounds checking
//     cout << "\nAccessing out-of-bounds elements:" << endl;
//     cout << "Attempt to access v10[-1]" << endl;
//     double value = v10[-1];  // Should print error message

//     cout << "Attempt to access v10(0)" << endl;
//     value = v10(0);  // Should print error message

//     return 0;
// }

