#ifndef MATRIX_H
#define MATRIX_H

#include <map>
#include <array> // per std::pair
#include <iostream>
#include <vector>
// Libraries for read function
#include <istream>      // per std::istream
#include <string>       // per std::string
#include <sstream>      // per std::istringstream
#include <fstream>      // per std::ifstream

// Define the algebra namespace
namespace algebra {

    // Define a comparison structure for sorting the storage
    template <std::size_t StorageOrder>
    struct compare {
        // Define the comparison operator
        bool operator()(const std::array<std::size_t, 2>& lhs, const std::array<std::size_t, 2>& rhs) const;
    };

    // Define the Matrix class
    template <typename T, std::size_t StorageOrder>  // StorageOrder = 0 for row-major, 1 for column-major
    class Matrix {

        // Define necessary data types
        typedef std::map<std::array<std::size_t, 2>, T, compare<StorageOrder>> MatrixMap;
        MatrixMap data;
        std::size_t rows=0;
        std::size_t cols=0;

        std::vector<T> Values;
        std::vector<std::size_t> Inner_dir;
        std::vector<std::size_t> Outer_dir;
        bool Compressed=false;

    public:
        // Define constructors
        Matrix() = default;
        Matrix(const MatrixMap& inputData) : data(inputData) {
            for (const auto& [key, value] : inputData) {
            rows = std::max(rows, key[0] + 1);
            cols = std::max(cols, key[1] + 1);
        }
        };
        Matrix(const Matrix& other) : data(other.data), rows(other.rows), cols(other.cols) {};
        Matrix(std::size_t R, std::size_t C) : rows(R), cols(C) {};

        // Define methods
        void resize(std::size_t R, std::size_t C);

        T& operator()(const std::array<std::size_t, 2>& index);
        T operator()(const std::array<std::size_t, 2>& index) const;
        // Easier access for the call operator
        T& operator()(const std::size_t& row, const std::size_t& col){return this->operator()({row, col});}
        T operator()(const std::size_t& row, const std::size_t& col) const {return this->operator()({row, col});}
        // Call for a row or a column based on the StorageOrder
        std::vector<T> operator () (std::size_t index);

        void compress();
        void uncompress();

        void read(const std::string& filename);
        void print() const;

        // Getters for the dimensions
        std::size_t getRows() const { return rows; }
        std::size_t getCols() const { return cols; }

        // Method to check if the matrix is compressed
        bool is_compressed() const { return Compressed; }

        // Friend functions for multiplication
        friend std::vector<T> operator* <> (const Matrix<T,StorageOrder>& lhs, const std::vector<T>& rhs);

        // Fix the StorageOrder of the result matrix as row-wise (for column-wise the functions are similar)
        template <typename U, std::size_t SO1, std::size_t SO2>
        friend Matrix<U,0> operator*(const Matrix<U, SO1>& lhs, const Matrix<U, SO2>& rhs);
    };
}

#include "Matrix_impl.hpp"

#endif // MATRIX_H