#include "Matrix.hpp"
#include <iostream>

#include <chrono>


int main(){
    /*  -----  READING FROM MM FILE  -----  */
    algebra::Matrix<double, 0> Matrix;
    std::string filename = "lnsp_131.mtx";
    Matrix.read(filename);

    /*  -----  MATRIX-VECTOR MULTIPLICATION  -----  */

    std::cout << "\nMatrix-Vector Multiplication:" << std::endl;
    std::vector<double> v(Matrix.getCols(),5);

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> result = Matrix * v;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << "Time to multiply: " << diff.count() << " s"<< std::endl;

    Matrix.compress();
    start = std::chrono::high_resolution_clock::now();
    result = Matrix * v;
    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Time to multiply after compression: " << diff.count() << " s\n"<< std::endl;

    /*  -----  MATRIX-MATRIX MULTIPLICATION  -----  */
    Matrix.uncompress();
    algebra::Matrix<double, 0> Matrix2(Matrix);
    Matrix2.resize(Matrix.getCols(),3);

    std::cout << "\nMatrix-Matrix Multiplication:" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    algebra::Matrix<double, 0> result2 = Matrix * Matrix2;
    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Time to multiply: " << diff.count() << " s"<< std::endl;

    Matrix.compress();
    Matrix2.compress();
    start = std::chrono::high_resolution_clock::now();
    result2 = Matrix * Matrix2;
    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Time to multiply after compression: " << diff.count() << " s"<< std::endl;
    
    return 0;
}