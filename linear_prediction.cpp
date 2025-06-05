#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>

#include "matrix.h"
#include "vector.h"

using namespace std;

void printVector(const Vector& v, string name); // function outside class vector, need to declare before using
void printMatrix(const Matrix& m, string name){
    cout << name << ":\n";
    for(int i = 0; i < m.getNumRows(); i++){
        for(int j = 0; j < m.getNumCols(); j++){
            cout << m.getData(i, j) << " ";
        }
        cout << endl;
    }
};

// Fix the shuffle function to use modern C++ random facilities correctly
void shuffle(vector<string>& lines) {
    // Create a random device and generator
    random_device rd;
    mt19937 g(rd());  // Mersenne Twister generator
    
    //Standard library shuffle algorithm
    shuffle(lines.begin(), lines.end(), g);
}

void loadInput(ifstream& datafile, Matrix& m){
    string line;

    // datafile.clear();
    // datafile.seekg(0);

    int rowIndex = 0;
    while(getline(datafile, line)){
        stringstream ss(line);
        string value;
        int index = 0;
        int colIndex = 0;

        while(getline(ss, value, ',')){
            if(index >= 2 && index <= 8){
                double d_value = stod(value);
                m(rowIndex + 1, colIndex + 1) = d_value; 
                colIndex++;
            }
            index++;
        }
        rowIndex++;

        if(rowIndex >= m.getNumRows()){
            break;
        }
    }
};
void loadOutput(ifstream& datafile, Matrix& m){
    string line;

    // datafile.clear();
    // datafile.seekg(0);

    int rowIndex = 0;
    while(getline(datafile, line)){
        stringstream ss(line);
        string value;
        int index = 0;

        while(getline(ss, value, ',')){
            if(index == 9){
                double d_value = stod(value);
                m(rowIndex + 1, 1) = d_value;
            }
            index++;
        }
        rowIndex++;

        if(rowIndex >= m.getNumRows()){
            break;
        }
    }
};

void normalizeMatrix(Matrix& m) {
    // For each column
    for(int j = 0; j < m.getNumCols(); j++) {
        // Find mean and std for column
        double sum = 0.0, sqSum = 0.0;
        for(int i = 0; i < m.getNumRows(); i++) {
            sum += m.getData(i, j);
            sqSum += m.getData(i, j) * m.getData(i, j);
        }
        double mean = sum / m.getNumRows();
        double std = sqrt((sqSum / m.getNumRows()) - (mean * mean));
        
        // Skip if std is too small (constant column)
        if(std < 1e-10) continue;
        
        // Normalize: (x - mean) / std
        for(int i = 0; i < m.getNumRows(); i++) {
            m(i+1, j+1) = (m.getData(i, j) - mean) / std;
        }
    }
}

int main(){
    ifstream Dataset("machine.data");

    int counter = 0;

    string line;
    vector<string> lines;
    while(std::getline(Dataset, line)){
        counter += 1;
        lines.push_back(line);
    }
    Dataset.close();

    shuffle(lines);
    ofstream outfile("shuffled.txt"); // shuffle machine.data and put to shiffled.txt
    for(const auto& l: lines){
        outfile << l << endl;
    }
    outfile.close();

    int trainSet = counter * 0.8;
    int testSet = counter - trainSet;
    
    cout << "Number of data(lines): " << counter << endl;
    cout << "Training set size: " << trainSet << endl;
    cout << "Testing set size: " << testSet << endl;

    ifstream shuffle_dataset("shuffled.txt");

    // Create matrices for training 
    Matrix trainInput(trainSet, 7);
    Matrix trainTarget(trainSet, 1); 

    loadInput(shuffle_dataset, trainInput);
    loadOutput(shuffle_dataset, trainTarget);

    //printMatrix(trainInput, "Training Input Matrix");
    //printMatrix(trainTarget, "Training Target Matrix");

    shuffle_dataset.clear();
    shuffle_dataset.seekg(0);

    // Create matrices for testing
    Matrix testInput(testSet, 7);
    Matrix testTarget(testSet, 1);

    string linee;
    for(int i = 0; i < trainSet; i++){
        getline(shuffle_dataset, linee);
    }

    loadInput(shuffle_dataset, testInput);
    loadOutput(shuffle_dataset, testTarget);

    shuffle_dataset.close();

    normalizeMatrix(trainInput);
    normalizeMatrix(testInput);

    // printMatrix(testInput, "Test input");
    // printMatrix(testTarget, "Test target");

    // Calculate theta: θ = (Xt*X)^(-1) * (Xt*Y)
    Matrix Xt = trainInput.tranpose();
    Matrix XtX = Xt * trainInput;

    // Apply regularization to improve numerical stability
    double lambda = 10.0;  // Regularization parameter (10.0 currently have the best RMSE value)
    for(int i = 1; i <= XtX.getNumRows(); i++) {
        XtX(i, i) += lambda;
    }
    cout << "Applied regularization with lambda = " << lambda << endl;

    Matrix XtX_inverse = XtX.pseudo_inverse();
    printMatrix(XtX_inverse, "XtX:");
    Matrix Xty = Xt * trainTarget;
    Matrix theta = XtX_inverse * Xty; // return a 7 rows matrix (6 coefficients and 1 bias) 
    //printMatrix(theta, "Model paremeters:");

    Matrix train_prediction = trainInput * theta;
    //printMatrix(train_prediction, "Train prediction:");

    double train_RMSE = 0.0;
    for(int i =0; i < trainTarget.getNumRows(); i++){
        double sum = trainTarget.getData(i, 0) - train_prediction.getData(i, 0);
        train_RMSE += sum * sum;
    }

    train_RMSE = sqrt(train_RMSE / trainTarget.getNumRows());
    cout << "Train RMSE: " << train_RMSE << endl;

    Matrix test_prediction = testInput * theta;
    //printMatrix(test_prediction, "Test prediction");
    
    double test_RMSE = 0.0;
    for(int i = 0; i < testTarget.getNumRows(); i++){
        double x = testTarget.getData(i, 0) - test_prediction.getData(i, 0);
        test_RMSE += x * x;
    }

    test_RMSE = sqrt(test_RMSE/ testTarget.getNumRows());
    cout << "Test RMSE: " << test_RMSE << endl;

    
    cout << "\nLinear Regression Model:" << endl;
    cout << "PRP = ";
    for (int i = 0; i < theta.getNumRows(); i++) {
        if (i == 0) {
            cout << theta.getData(i, 0) << " * MYCT ";
        } else if (i == 1) {
            cout << (theta.getData(i, 0) >= 0 ? "+ " : "- ") << abs(theta.getData(i, 0)) << " * MMIN ";
        } else if (i == 2) {
            cout << (theta.getData(i, 0) >= 0 ? "+ " : "- ") << abs(theta.getData(i, 0)) << " * MMAX ";
        } else if (i == 3) {
            cout << (theta.getData(i, 0) >= 0 ? "+ " : "- ") << abs(theta.getData(i, 0)) << " * CACH ";
        } else if (i == 4) {
            cout << (theta.getData(i, 0) >= 0 ? "+ " : "- ") << abs(theta.getData(i, 0)) << " * CHMIN ";
        } else if (i == 5) {
            cout << (theta.getData(i, 0) >= 0 ? "+ " : "- ") << abs(theta.getData(i, 0)) << " * CHMAX ";
        } else if (i == 6) {
            cout << (theta.getData(i, 0) >= 0 ? "+ " : "- ") << abs(theta.getData(i, 0));
        }
    }
    cout << endl;

    return 0;
}

// g++ -o test linear_prediction.cpp "Part A/matrix.cpp" "Part A/vector.cpp" -I"Part A"
// ./test