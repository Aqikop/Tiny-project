#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "matrix.h"
#include "vector.h"

using namespace std;

void printVector(const Vector& v, string name); // function outside class vector, need to declare before using

int main(){
    ifstream Dataset("machine.data");

    int count = 0;

    int arr[7] = {0};

    string line;
    while(std::getline(Dataset, line)){
        count += 1;
        
        std::stringstream a(line);
        string value;
        int index = 0;
        
        while(std::getline(a, value, ',')){
            if(index >=2 && index <= 8){
                int intValue = std::stoi(value); // Convert str to int
                arr[index - 2] = intValue;
            }
            index++;
        }
        break;
    }
    
    // cout << count << endl;

    // int trainSet = count * 0.8;
    // int testSet = count - trainSet;

    Vector v(7);
    for(int i = 0; i < v.get_size(); i++){
        v[i] = arr[i];
    }
    printVector(v, "First line");

    Dataset.close();
    return 0;
}

// Compile with: g++ -std=c++17 -Wall -I"Part A" test.cpp "Part A/vector.cpp" "Part A/matrix.cpp" -o test
// Run with: ./test
// Delete object file with: rm test
